# schur.jl -- the fixed-effect solver of the v1.0.0 fast path (solver = :schur).
#
# Borrowed from fereg's fedemean (Schmieder 2026, NBER WP 36017, p.22-24): eliminate the fixed
# effect with the most groups EXACTLY, then solve the remaining, much smaller system (the Schur
# complement S = N2 - C' N1^{-1} C, a graph-Laplacian-like matrix) by Jacobi-preconditioned
# conjugate gradients, all data columns advancing in lockstep through one pass over the
# observations. S is never formed: S*x = D2' W M1 D2 x, where M1 is the weighted within-group
# demeaning for the eliminated FE. With one FE this is the exact weighted within-group demeaning.
#
# Stopping rule (per column): ||A'e||_2 <= tol * sqrt(F) * ||sqrt(W) e||_2, with
# A = sqrt(W) D diag(1/sqrt(group weight)) the column-normalised design that FixedEffects.jl's LSMR
# works on. ||A||_2 <= sqrt(F) exactly (each FE block of A has orthonormal columns), so this is never
# looser than LSMR's own test2 (normAr/(normA*normr) <= atol, with normA its running estimate, which
# grows with the iteration count). The FE-1 part of A'e is zero by construction (a1 is solved
# exactly), so only the remaining levels enter. A column that passes on the recursive CG residual
# must pass again on a freshly recomputed residual (a certificate pass, as in fereg) before it is
# accepted; otherwise CG restarts from the current iterate.
#
# Disconnected graphs: S is singular but the system is consistent; CG on the range is fine, and
# levels whose Jacobi diagonal vanishes (absorbed by FE 1) get a zero preconditioner entry.
#
# Origin: program/sandbox/julia/fast_demean_20261119/scripts/IVPPMLFast/src/schur.jl (2026-11-19),
# unit-tested there against LSMR at 1e-14 (residuals agree to 2e-16 .. 2e-12 relative).

mutable struct SchurFESolver
    N::Int
    F::Int                       # number of FEs
    big::Int                     # index of the eliminated FE (most levels)
    G::Int                       # its number of groups
    perm::Vector{Int}            # sorted position k -> original observation perm[k]
    isid::Bool                   # perm is the identity
    gptr::Vector{Int}            # group g occupies sorted positions gptr[g]:gptr[g+1]-1
    maxg::Int                    # largest group size
    R::Int                       # number of remaining FEs
    lev::Matrix{Int32}           # R x N remaining-level index in 1..m (sorted order)
    m::Int                       # number of remaining levels
    w::Vector{Float64}           # weights in sorted order
    n1inv::Vector{Float64}       # 1 / (group weight sum) for the eliminated FE
    d2::Vector{Float64}          # weight sum per remaining level
    d2inv::Vector{Float64}
    pinv::Vector{Float64}        # inverse Jacobi diagonal of S (0 for levels absorbed by FE 1)
    sqrtF::Float64
    nc::Int                      # buffer width (number of columns)
    V::Matrix{Float64}           # nc x N  input data (sorted order)
    E::Matrix{Float64}           # nc x N  residuals (sorted order)
    X::Matrix{Float64}           # nc x m  CG iterate
    Rr::Matrix{Float64}          # nc x m  CG residual
    P::Matrix{Float64}           # nc x m  search direction
    Q::Matrix{Float64}           # nc x m  S * P  (also: recomputed gradient)
    B::Matrix{Float64}           # nc x m  right-hand side
    U::Matrix{Float64}           # nc x maxg  per-group scratch
    acc::Vector{Float64}         # length m scratch
    sbuf::Vector{Float64}        # length nc scratch
    # counters (cumulative)
    passes::Int                  # passes over the data (operator applications + residual passes)
    restarts::Int                # CG restarts forced by a failed certificate
end

function SchurFESolver(fes::AbstractVector, w::AbstractVector{<:Real})
    F = length(fes)
    N = length(fes[1].refs)
    for fe in fes
        (fe.interaction isa UnitWeights || all(==(1), fe.interaction)) ||
            error("SchurFESolver: interacted (slope) fixed effects are not supported")
    end
    big = argmax([fe.n for fe in fes])           # ties: the first
    refs1 = fes[big].refs
    G = fes[big].n
    # counting sort of observations by the eliminated FE (stable)
    cnt = zeros(Int, G + 1)
    @inbounds for i in 1:N
        cnt[refs1[i] + 1] += 1
    end
    gptr = Vector{Int}(undef, G + 1)
    gptr[1] = 1
    @inbounds for g in 1:G
        gptr[g + 1] = gptr[g] + cnt[g + 1]
    end
    pos = gptr[1:G]
    perm = Vector{Int}(undef, N)
    @inbounds for i in 1:N
        g = refs1[i]
        perm[pos[g]] = i
        pos[g] += 1
    end
    isid = all(k -> perm[k] == k, 1:N)
    maxg = maximum(gptr[g + 1] - gptr[g] for g in 1:G)
    rest = [f for f in 1:F if f != big]
    R = length(rest)
    lev = Matrix{Int32}(undef, R, N)
    off = 0
    for (r, f) in enumerate(rest)
        rf = fes[f].refs
        @inbounds for k in 1:N
            lev[r, k] = Int32(off + rf[perm[k]])
        end
        off += fes[f].n
    end
    m = off
    S = SchurFESolver(N, F, big, G, perm, isid, gptr, maxg, R, lev, m,
                      zeros(N), zeros(G), zeros(m), zeros(m), zeros(m), sqrt(Float64(F)),
                      0, zeros(0, 0), zeros(0, 0), zeros(0, 0), zeros(0, 0), zeros(0, 0),
                      zeros(0, 0), zeros(0, 0), zeros(0, 0), zeros(m), zeros(0), 0, 0)
    schur_update_weights!(S, w)
    return S
end

function _schur_buffers!(S::SchurFESolver, nc::Int)
    S.nc == nc && return S
    S.nc = nc
    S.V = zeros(nc, S.N); S.E = zeros(nc, S.N)
    S.X = zeros(nc, S.m); S.Rr = zeros(nc, S.m); S.P = zeros(nc, S.m)
    S.Q = zeros(nc, S.m); S.B = zeros(nc, S.m)
    S.U = zeros(nc, S.maxg); S.sbuf = zeros(nc)
    return S
end

# Refill the weight-dependent pieces only (the structure is built once).
function schur_update_weights!(S::SchurFESolver, w::AbstractVector{<:Real})
    N, G, R, m = S.N, S.G, S.R, S.m
    perm, gptr, lev = S.perm, S.gptr, S.lev
    sw, n1inv, d2, acc, pinv = S.w, S.n1inv, S.d2, S.acc, S.pinv
    @inbounds for k in 1:N
        sw[k] = w[perm[k]]
    end
    fill!(d2, 0.0)
    @inbounds for g in 1:G
        s = 0.0
        for k in gptr[g]:(gptr[g + 1] - 1)
            wk = sw[k]
            s += wk
            for r in 1:R
                d2[lev[r, k]] += wk
            end
        end
        n1inv[g] = s > 0 ? 1.0 / s : 0.0
    end
    # diag(S)_l = d2_l - sum_g C_{g,l}^2 / n1_g ,  C_{g,l} = sum_{i in g, lev(i)=l} w_i
    copyto!(pinv, d2)                 # pinv holds diag(S) for now
    fill!(acc, 0.0)
    @inbounds for g in 1:G
        lo, hi = gptr[g], gptr[g + 1] - 1
        for k in lo:hi, r in 1:R
            acc[lev[r, k]] += sw[k]
        end
        ni = n1inv[g]
        for k in lo:hi, r in 1:R
            l = lev[r, k]
            a = acc[l]
            if a != 0.0
                pinv[l] -= a * a * ni
                acc[l] = 0.0
            end
        end
    end
    @inbounds for l in 1:m
        dS = pinv[l]
        pinv[l] = dS > 1e-12 * d2[l] ? 1.0 / dS : 0.0
        S.d2inv[l] = d2[l] > 0 ? 1.0 / d2[l] : 0.0
    end
    return S
end

# Q = S * P  for all nc columns in one pass over the data. The kernels are specialised on the
# number of columns and of remaining FEs (compile-time loop bounds, so the inner loops unroll).
function _schur_apply!(S::SchurFESolver, Q::Matrix{Float64}, P::Matrix{Float64})
    _schur_apply_k!(S, Q, P, Val(S.nc), Val(S.R))
    S.passes += 1
    return Q
end

function _schur_apply_k!(S::SchurFESolver, Q::Matrix{Float64}, P::Matrix{Float64},
                         ::Val{nc}, ::Val{R}) where {nc, R}
    U, sb, w, lev, gptr = S.U, S.sbuf, S.w, S.lev, S.gptr
    fill!(Q, 0.0)
    @inbounds for g in 1:S.G
        lo, hi = gptr[g], gptr[g + 1] - 1
        for c in 1:nc
            sb[c] = 0.0
        end
        for k in lo:hi
            j = k - lo + 1
            for c in 1:nc
                U[c, j] = 0.0
            end
            for r in 1:R
                l = lev[r, k]
                for c in 1:nc
                    U[c, j] += P[c, l]
                end
            end
            wk = w[k]
            for c in 1:nc
                sb[c] += wk * U[c, j]
            end
        end
        ni = S.n1inv[g]
        for c in 1:nc
            sb[c] *= ni
        end
        for k in lo:hi
            j = k - lo + 1
            wk = w[k]
            for r in 1:R
                l = lev[r, k]
                for c in 1:nc
                    Q[c, l] += wk * (U[c, j] - sb[c])
                end
            end
        end
    end
    return Q
end

# E = M1 (V - D2 X);  Gr = D2' W E;  en2[c] = ||sqrt(W) E[c,:]||^2.   (with X = 0: Gr = rhs, en2 = c0)
function _schur_residual!(S::SchurFESolver, Gr::Matrix{Float64}, en2::Vector{Float64},
                          X::Matrix{Float64}, writeE::Bool)
    _schur_residual_k!(S, Gr, en2, X, writeE, Val(S.nc), Val(S.R))
    S.passes += 1
    return Gr
end

function _schur_residual_k!(S::SchurFESolver, Gr::Matrix{Float64}, en2::Vector{Float64},
                            X::Matrix{Float64}, writeE::Bool, ::Val{nc}, ::Val{R}) where {nc, R}
    U, sb, w, lev, gptr, V, E = S.U, S.sbuf, S.w, S.lev, S.gptr, S.V, S.E
    fill!(Gr, 0.0)
    fill!(en2, 0.0)
    @inbounds for g in 1:S.G
        lo, hi = gptr[g], gptr[g + 1] - 1
        for c in 1:nc
            sb[c] = 0.0
        end
        for k in lo:hi
            j = k - lo + 1
            for c in 1:nc
                U[c, j] = V[c, k]
            end
            for r in 1:R
                l = lev[r, k]
                for c in 1:nc
                    U[c, j] -= X[c, l]
                end
            end
            wk = w[k]
            for c in 1:nc
                sb[c] += wk * U[c, j]
            end
        end
        ni = S.n1inv[g]
        for c in 1:nc
            sb[c] *= ni
        end
        for k in lo:hi
            j = k - lo + 1
            wk = w[k]
            for c in 1:nc
                e = U[c, j] - sb[c]
                U[c, j] = wk * e
                en2[c] += wk * e * e
                writeE && (E[c, k] = e)
            end
            for r in 1:R
                l = lev[r, k]
                for c in 1:nc
                    Gr[c, l] += U[c, j]
                end
            end
        end
    end
    return Gr
end

# Partial the fixed effects out of every column of `data` (N x nc), in place, under the current
# weights. Returns (data, iterations per column, converged per column).
function schur_solve_residuals!(data::AbstractMatrix{Float64}, S::SchurFESolver;
                                tol::Real = 1e-8, maxiter::Integer = 10000)
    N, nc = size(data)
    N == S.N || error("schur_solve_residuals!: row mismatch")
    _schur_buffers!(S, nc)
    V, E, X, Rr, P, Q, B, d2inv, pinv, perm, m = S.V, S.E, S.X, S.Rr, S.P, S.Q, S.B, S.d2inv, S.pinv, S.perm, S.m
    @inbounds if S.isid
        for c in 1:nc, k in 1:N
            V[c, k] = data[k, c]
        end
    else
        for c in 1:nc, k in 1:N
            V[c, k] = data[perm[k], c]
        end
    end
    thr2 = (Float64(tol) * S.sqrtF)^2
    fill!(X, 0.0)
    c0 = zeros(nc)
    _schur_residual!(S, B, c0, X, false)          # B = rhs, c0 = ||sqrt(W) M1 v||^2
    copyto!(Rr, B)
    iters = zeros(Int, nc)
    conv = falses(nc)
    active = trues(nc)
    rz = zeros(nc)
    en2 = zeros(nc)
    gs(c, A) = (s = 0.0; @inbounds for l in 1:m; s += A[c, l]^2 * d2inv[l]; end; s)
    # columns already converged at x = 0
    for c in 1:nc
        if gs(c, Rr) <= thr2 * c0[c]
            active[c] = false
        end
    end
    total = 0
    restarts = 0
    while true
        # (re)start CG for the active columns from the current X with residual Rr
        @inbounds for c in 1:nc
            if active[c]
                s = 0.0
                for l in 1:m
                    z = pinv[l] * Rr[c, l]
                    P[c, l] = z
                    s += Rr[c, l] * z
                end
                rz[c] = s
                rz[c] > 0 || (active[c] = false)
            else
                for l in 1:m
                    P[c, l] = 0.0
                end
            end
        end
        while any(active) && total < maxiter
            _schur_apply!(S, Q, P)
            total += 1
            @inbounds for c in 1:nc
                active[c] || continue
                pq = 0.0
                for l in 1:m
                    pq += P[c, l] * Q[c, l]
                end
                if !(pq > 0)
                    active[c] = false
                    for l in 1:m; P[c, l] = 0.0; end
                    continue
                end
                α = rz[c] / pq
                xb = 0.0
                for l in 1:m
                    X[c, l] += α * P[c, l]
                    Rr[c, l] -= α * Q[c, l]
                    xb += X[c, l] * (B[c, l] + Rr[c, l])
                end
                iters[c] += 1
                e2 = max(c0[c] - xb, 0.0)            # ||sqrt(W) e||^2 = c0 - x'(b + r)
                if gs(c, Rr) <= thr2 * e2
                    active[c] = false
                    for l in 1:m; P[c, l] = 0.0; end
                    continue
                end
                s = 0.0
                for l in 1:m
                    s += Rr[c, l] * pinv[l] * Rr[c, l]
                end
                β = s / rz[c]
                rz[c] = s
                for l in 1:m
                    P[c, l] = pinv[l] * Rr[c, l] + β * P[c, l]
                end
            end
        end
        # certificate: recompute the residual and the gradient from scratch
        _schur_residual!(S, Q, en2, X, true)
        again = false
        @inbounds for c in 1:nc
            conv[c] && continue
            if gs(c, Q) <= thr2 * en2[c]
                conv[c] = true
            elseif total < maxiter && restarts < 50
                for l in 1:m
                    Rr[c, l] = Q[c, l]
                end
                active[c] = true
                again = true
            end
        end
        again || break
        restarts += 1
    end
    S.restarts += restarts
    @inbounds if S.isid
        for c in 1:nc, k in 1:N
            data[k, c] = E[c, k]
        end
    else
        for c in 1:nc, k in 1:N
            data[perm[k], c] = E[c, k]
        end
    end
    return data, iters, collect(conv)
end
