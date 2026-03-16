# IVPPMLFixedEffectModels.jl — IV-PPML with High-Dimensional Fixed Effects
#
# Path 3 architecture: single-call estimation with solver reuse,
# mirroring GLFixedEffectModels.jl's nlreg() pattern.
#
# Key difference from ivppmlhdfe.jl:
#   - FE solver created ONCE, weights updated via update_weights!()
#   - DataFrame-based input (like reghdfejl)
#   - Result struct with standard accessors (coef, vcov, nobs, etc.)
#
# Algorithm: IRLS-IV (Mullahy 1997 GMM via iteratively reweighted 2SLS)
# FE absorption via FixedEffects.jl (iterative demeaning / FWL)
#
# Authors: Ohyun Kwon
# Date: 2026-03-14

module IVPPMLFixedEffectModels

using LinearAlgebra
using FixedEffects
using FixedEffects: update_weights!
using StatsBase: Weights
using Printf
using PrecompileTools

export ivppml_reg, IVPPMLResult, coef, vcov, nobs, coefnames

# ============================================================================
# Result container with standard accessors
# ============================================================================

struct IVPPMLResult
    coef_vec::Vector{Float64}
    vcov_mat::Matrix{Float64}
    coefnames_vec::Vector{String}
    N::Int
    N_full::Int
    converged::Bool
    iterations::Int
    deviance::Float64
    esample::BitVector
end

# Standard accessors (matching FixedEffectModels / GLFixedEffectModels API)
coef(r::IVPPMLResult)      = r.coef_vec
vcov(r::IVPPMLResult)      = r.vcov_mat
nobs(r::IVPPMLResult)      = r.N
coefnames(r::IVPPMLResult) = r.coefnames_vec

# ============================================================================
# Weighted 2SLS on demeaned data (no constant)
# ============================================================================

function _ivppml_2sls!(b::Vector{Float64},
                       resid::Vector{Float64},
                       Xhat::Matrix{Float64},
                       z_dm::AbstractVector{Float64},
                       X_dm::AbstractMatrix{Float64},
                       Z_dm::AbstractMatrix{Float64},
                       w::AbstractVector{Float64})
    # First stage:  Pi = (Z'WZ)^{-1} Z'WX
    ZwZ = Z_dm' * (w .* Z_dm)
    ZwX = Z_dm' * (w .* X_dm)
    Pi  = ZwZ \ ZwX
    Xhat .= Z_dm * Pi

    # Second stage: b = (Xhat'WX)^{-1} Xhat'Wy
    XhwX = Xhat' * (w .* X_dm)
    Xhwy = Xhat' * (w .* z_dm)
    b   .= XhwX \ Xhwy

    resid .= z_dm .- X_dm * b
    return nothing
end

# ============================================================================
# Robust sandwich VCE  (HC1: N/(N-K) small-sample correction)
# ============================================================================

function _ivppml_robust_vce(Xhat::Matrix{Float64},
                            X_dm::Matrix{Float64},
                            w::Vector{Float64},
                            resid::Vector{Float64},
                            N::Int)
    K = size(X_dm, 2)
    bread  = inv(Symmetric(Xhat' * (w .* X_dm)))
    scores = Xhat .* (w .* resid)
    meat   = scores' * scores
    V = (N / (N - K)) .* bread * meat * bread
    return V
end

# ============================================================================
# Cluster-robust sandwich VCE  (Arellano G/(G-1) correction)
# ============================================================================

function _ivppml_cluster_vce(Xhat::Matrix{Float64},
                             X_dm::Matrix{Float64},
                             w::Vector{Float64},
                             resid::Vector{Float64},
                             clust_id::AbstractVector)
    K = size(X_dm, 2)
    bread  = inv(Symmetric(Xhat' * (w .* X_dm)))
    scores = Xhat .* (w .* resid)

    perm = sortperm(clust_id)
    sorted_id = clust_id[perm]
    sorted_scores = scores[perm, :]

    G = 0
    meat = zeros(K, K)
    i = 1
    while i <= length(sorted_id)
        j = i
        while j <= length(sorted_id) && sorted_id[j] == sorted_id[i]
            j += 1
        end
        sg = vec(sum(view(sorted_scores, i:j-1, :); dims=1))
        meat .+= sg * sg'
        G += 1
        i = j
    end

    V = (G / (G - 1)) .* bread * meat * bread
    return V, G
end

# ============================================================================
# Main estimation: ivppml_reg
#
# Single-call entry point mirroring nlreg() from GLFixedEffectModels.jl.
# Accepts a DataFrame and column symbols.
# ============================================================================

function ivppml_reg(df::Any,          # AbstractDataFrame
                    depvar::Symbol,
                    exog::Vector{Symbol},
                    endog::Vector{Symbol},
                    instruments::Vector{Symbol},
                    fe_syms::Vector{Symbol};
                    cluster::Union{Symbol, Nothing} = nothing,
                    weights::Union{Symbol, Nothing} = nothing,
                    tol::Float64  = 1e-8,
                    maxiter::Int  = 1000,
                    verbose::Int  = 0,
                    method::Symbol = :cpu)

    # ---- Extract data from DataFrame ----
    N_full = size(df, 1)
    y = Float64.(df[!, depvar])
    N = length(y)

    K_exog  = length(exog)
    K_endo  = length(endog)
    K       = K_exog + K_endo
    n_inst  = length(instruments)
    L       = K_exog + n_inst         # total instruments

    # Build X = [exog endog]
    X = Matrix{Float64}(undef, N, K)
    for (j, s) in enumerate(exog)
        X[:, j] .= Float64.(df[!, s])
    end
    for (j, s) in enumerate(endog)
        X[:, K_exog + j] .= Float64.(df[!, s])
    end

    # Build excluded instruments
    Z_excl = Matrix{Float64}(undef, N, n_inst)
    for (j, s) in enumerate(instruments)
        Z_excl[:, j] .= Float64.(df[!, s])
    end

    # Build FE objects
    fes = [FixedEffect(df[!, s]) for s in fe_syms]

    # User weights
    w_user = weights === nothing ? ones(Float64, N) : Float64.(df[!, weights])

    # Cluster IDs
    has_cluster = cluster !== nothing
    clust_id = has_cluster ? df[!, cluster] : nothing

    # Coefficient names
    cnames = vcat(String.(vcat(exog, endog)), "_cons")

    # ---- Initialise mu ----
    sw = sum(w_user)
    mean_y = dot(w_user, y) / sw
    mu  = @. 0.5 * (y + mean_y)
    mu  = max.(mu, 1e-4)
    eta = log.(mu)

    # ---- Build FE solver ONCE ----
    irls_w = w_user .* mu
    wts    = Weights(irls_w)
    feM    = AbstractFixedEffectSolver{Float64}(fes, wts, Val{method})

    # ---- IRLS state ----
    converged = false
    deviance  = Inf
    eps_val   = Inf
    ok        = 0
    iter      = 0

    # Adaptive inner tolerance (following ppmlhdfe)
    start_inner_tol  = 1e-4
    target_inner_tol = max(1e-12, 0.1 * tol)
    inner_tol = max(start_inner_tol, tol)
    alt_tol   = start_inner_tol

    # Step-halving
    max_step_halving   = 2
    step_memory        = 0.9
    num_step_halving   = 0

    # Pre-allocate
    b     = zeros(K)
    resid = zeros(N)
    Xhat  = zeros(N, K)
    old_eta = similar(eta)

    ncols = 1 + K + n_inst
    data  = Matrix{Float64}(undef, N, ncols)
    X_dm  = Matrix{Float64}(undef, N, K)
    Z_dm  = Matrix{Float64}(undef, N, L)

    if verbose > -1
        @printf("IRLS-IV iterations (N = %d, K = %d, L = %d)\n", N, K, L)
        println("-"^60)
    end

    for it in 1:maxiter
        iter = it

        # (a) Working dependent variable: z = eta - 1 + y / mu
        z = eta .- 1.0 .+ y ./ mu

        # (b) IRLS weights: w = w_user * mu
        @. irls_w = w_user * mu

        # (c) Update solver weights IN-PLACE (not recreate!)
        update_weights!(feM, Weights(irls_w))

        # (d) Demean all variables at once
        data[:, 1]            .= z
        data[:, 2:K+1]        .= X
        data[:, K+2:end]      .= Z_excl
        solve_residuals!(data, feM; maxiter=10000, tol=inner_tol, progress_bar=false)

        z_dm = view(data, :, 1)
        X_dm .= view(data, :, 2:K+1)

        # Build Z_dm = [exog_dm, inst_excl_dm]
        if K_exog > 0 && n_inst > 0
            Z_dm .= hcat(view(data, :, 2:K_exog+1), view(data, :, K+2:K+n_inst+1))
        elseif K_exog > 0
            Z_dm .= view(data, :, 2:K_exog+1)
        else
            Z_dm .= view(data, :, K+2:K+n_inst+1)
        end

        # (e) Weighted 2SLS
        _ivppml_2sls!(b, resid, Xhat, z_dm, X_dm, Z_dm, irls_w)

        # (f) Update eta
        copyto!(old_eta, eta)
        @. eta = z - resid

        # (g) Update mu
        @. mu = exp(eta)
        @. mu = max(mu, 1e-10)

        # (h) Deviance
        old_deviance = deviance
        dev1 = dot(mu .- y, w_user)
        dev2 = 0.0
        @inbounds for i in 1:N
            if y[i] > 0.0
                dev2 += w_user[i] * y[i] * (log(y[i]) - eta[i])
            end
        end
        deviance = 2.0 * (dev1 + dev2)
        deviance = max(deviance, 0.0)

        # (i) Convergence + step-halving + adaptive tolerance
        is_step_halving = false
        if iter > 1
            delta_dev = old_deviance - deviance
            denom     = max(min(deviance, old_deviance), 0.1)
            eps_val   = abs(delta_dev) / denom

            if eps_val < tol
                if inner_tol <= 1.1 * target_inner_tol || length(fes) <= 1
                    ok += 1
                    ok >= 1 && (converged = true)
                end
            elseif delta_dev < 0 && num_step_halving < max_step_halving
                @. eta = step_memory * old_eta + (1 - step_memory) * eta
                @. mu  = exp(eta)
                @. mu  = max(mu, 1e-10)
                is_step_halving = true
                ok = 0
                num_step_halving += 1
            else
                ok = 0
                num_step_halving = 0
            end
        end

        if verbose > -1
            @printf("Iter %3d:  dev = %-11.5e", iter, deviance)
            iter > 1 && @printf("  eps = %-9.4e", eps_val)
            @printf("  tol = %5.0e", inner_tol)
            is_step_halving && print("  H")
            ok > 0          && print("  O")
            println()
        end

        if is_step_halving
            deviance = old_deviance
            continue
        end

        converged && break

        # Adaptive inner tolerance
        if iter > 1 && eps_val < inner_tol
            inner_tol = max(min(0.1 * inner_tol, alt_tol), target_inner_tol)
            alt_tol = 10.0^(-ceil(log10(1.0 / max(0.1 * eps_val, eps(1.0)))))
        end
    end

    if !converged && verbose > -1
        @printf("Warning: failed to converge in %d iterations (eps = %9.4e)\n",
                maxiter, eps_val)
    elseif verbose > -1
        @printf("Converged in %d iterations (tol = %9.4e)\n", iter, tol)
    end

    # ==================================================================
    # Final beta with constant, and VCE
    # ==================================================================

    @. irls_w = w_user * mu
    sw_f = sum(irls_w)

    mean_eta = dot(irls_w, eta) / sw_f
    mean_X   = vec(sum(irls_w .* X; dims=1)) ./ sw_f
    b_cons   = mean_eta - dot(mean_X, b)

    b_full  = vcat(b, b_cons)
    K_total = K + 1

    # Recompute Xhat for VCE
    ZwZ_f = Z_dm' * (irls_w .* Z_dm)
    Pi_f  = ZwZ_f \ (Z_dm' * (irls_w .* X_dm))
    Xhat  .= Z_dm * Pi_f

    # VCE
    N_clust = 0
    if has_cluster
        V_slope, N_clust = _ivppml_cluster_vce(Xhat, X_dm, irls_w, resid, clust_id)
    else
        V_slope = _ivppml_robust_vce(Xhat, X_dm, irls_w, resid, N)
    end

    # Expand to include constant
    V = zeros(K_total, K_total)
    V[1:K, 1:K] = V_slope
    V[K_total, K_total] = V_slope[1, 1]   # placeholder SE for _cons

    # esample (all true since data is pre-filtered by Stata)
    esample = trues(N_full)

    return IVPPMLResult(b_full, V, cnames, N, N_full,
                        converged, iter, deviance, esample)
end

# ============================================================================
# Precompile workload — runs during package precompilation so that Julia
# caches native code.  Eliminates JIT on first real call.
# ============================================================================

# Minimal wrapper to make a NamedTuple behave like a DataFrame for precompilation
# (avoids depending on DataFrames.jl just for precompilation)
struct _DFWrap{T}
    nt::T
end
Base.size(d::_DFWrap, i::Int) = i == 1 ? length(first(d.nt)) : error()
Base.getindex(d::_DFWrap, ::typeof(!), s::Symbol) = getfield(d.nt, s)

@compile_workload begin
    # Tiny dummy problem: N=20, K=1 exog, 1 endog, 1 instrument, 1 FE
    _N = 20
    _dfwrap = _DFWrap((;
        y    = exp.(randn(_N)) .+ 1.0,
        x1   = randn(_N),
        x2   = randn(_N),
        z1   = randn(_N),
        fe1  = repeat(1:4, 5),
        cid  = repeat(1:5, 4),
    ))

    # Run estimation (silent) — triggers compilation of all code paths
    ivppml_reg(_dfwrap, :y, Symbol[:x1], Symbol[:x2], Symbol[:z1], Symbol[:fe1];
               cluster = :cid, tol = 1e-4, maxiter = 5, verbose = -1)
end

end # module
