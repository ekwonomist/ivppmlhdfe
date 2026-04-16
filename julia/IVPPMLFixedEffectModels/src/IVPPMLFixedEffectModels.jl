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
# Algorithm: IRLS-IV targeting E[Z(y-mu)]=0 (Mullahy 1997) via iteratively reweighted 2SLS
# FE absorption via FixedEffects.jl (iterative demeaning / FWL)
#
# Authors: Ohyun Kwon
# Date: 2026-03-14
# Updated: 2026-04-12 — v0.9.2 feature parity with ivppmlhdfe.ado
#   * standardize option (X,Z column normalization + back-transform)
#   * mu-separation detection (accumulating mask, y=0 z-reset)
#   * guess(:mean|:simple) option for mu initialization
#   * actual_rank via collinearity check
#   * eta clipping at -10 on repeated step-halving
#   * PSD fix for multi-way clustering
#   * user-overridable inner_tol
# Updated: 2026-04-14 — v0.9.3 loud-failure parity with Stata ivppmlhdfe v0.9.3
#   * rc=9003 runaway divergence: Inf/NaN in mu/eta/irls_w and max|b|>1e6
#     now raise ErrorException instead of returning converged=false (no more
#     silent junk results on heavily-separated panels)
#   * rc=430 non-convergence: failing to hit tol within maxiter now raises
#     ErrorException; users must retry with higher maxiter or standardize=true
#   * fweight treated as pweight (documented limitation — Julia has no
#     weight-type API; users needing frequency weights should use Stata native)
# Updated: 2026-04-14 — v0.9.4 two-stage collinearity parity with Stata v0.9.4
#   * _select_not_collinear(): weighted pivoted-QR rank detection
#   * Pre-IRLS per-block collinearity pipeline on exog / endog / inst,
#     each partialled out against the FE structure (mirrors ppmlhdfe's
#     remove_collinears but applied three times — once per IV block)
#   * Replaces v0.9.2's silent-LU-pivot garbage on manually-expanded factor
#     dummies (T3.1 in 02_stress_julia.log v0.9.2 produced coefs 1e13 with
#     converged=false; v0.9.4 drops the redundant column and converges)
#   * Graceful re-check of the order condition after collinearity removal

module IVPPMLFixedEffectModels

using LinearAlgebra
using FixedEffects
using FixedEffects: update_weights!
using StatsBase: Weights
using SpecialFunctions: loggamma
using Printf
using PrecompileTools
using Random
# Extend StatsAPI generics rather than defining shadow functions locally.
# Without this import, `using IVPPMLFixedEffectModels, GLFixedEffectModels`
# causes name collisions and `coef(r)` throws UndefVarError in Main.
import StatsAPI: coef, vcov, nobs, coefnames, deviance, loglikelihood,
    nulldeviance, dof_residual, islinear, responsename,
    stderror, confint, coeftable
using StatsBase: CoefTable
import Distributions: Normal, cdf, quantile

export ivppml_reg, ppml_reg, IVPPMLResult
# Re-export the StatsAPI generics so that `using IVPPMLFixedEffectModels`
# still makes `coef(r)`, `vcov(r)`, etc. available in the caller's scope.
export coef, vcov, nobs, coefnames, deviance, loglikelihood,
    nulldeviance, dof_residual, stderror, confint, coeftable, responsename

# ============================================================================
# Result container with standard accessors
# ============================================================================

struct IVPPMLResult
    coef_vec::Vector{Float64}
    vcov_mat::Matrix{Float64}
    coefnames_vec::Vector{String}
    depvar::Symbol                  # response variable name (for responsename)
    N::Int
    N_full::Int
    converged::Bool
    iterations::Int
    deviance::Float64
    ll::Float64
    ll_0::Float64
    d_values::Vector{Float64}       # FE sum (eta - offset - X*b - b_cons)
    eta_values::Vector{Float64}     # linear predictor
    N_clust::Vector{Int}            # cluster counts per cluster variable
    esample::BitVector
    num_sep_mu::Int                 # observations dropped by mu-separation
    rank::Int                       # actual rank (= K minus collinear cols)
    df_absorb::Int                  # absorbed-FE degrees of freedom
    n_dropped_exog::Int             # columns dropped by two-stage collinearity
    n_dropped_endog::Int
    n_dropped_inst::Int
end

# StatsAPI-compatible accessors. Defensive `copy()` on mutable fields so
# callers cannot corrupt the result struct by mutating the returned view.
coef(r::IVPPMLResult)         = copy(r.coef_vec)
vcov(r::IVPPMLResult)         = copy(r.vcov_mat)
nobs(r::IVPPMLResult)         = r.N
coefnames(r::IVPPMLResult)    = copy(r.coefnames_vec)
deviance(r::IVPPMLResult)     = r.deviance
loglikelihood(r::IVPPMLResult)= r.ll
nulldeviance(r::IVPPMLResult) = -2 * r.ll_0
# dof_residual subtracts both regressor rank AND absorbed-FE df, matching the
# reghdfe/ppmlhdfe convention (e(df) = N - K - df_a).
dof_residual(r::IVPPMLResult) = max(r.N - r.rank - r.df_absorb, 0)
islinear(::IVPPMLResult)      = false
responsename(r::IVPPMLResult) = String(r.depvar)
stderror(r::IVPPMLResult)     = sqrt.(abs.(diag(r.vcov_mat)))

function confint(r::IVPPMLResult; level::Real=0.95)
    α = 1 - level
    z = quantile(Normal(0, 1), 1 - α/2)
    se = stderror(r)
    b = coef(r)
    return hcat(b .- z .* se, b .+ z .* se)
end

function coeftable(r::IVPPMLResult; level::Real=0.95)
    b  = coef(r)
    se = stderror(r)
    z  = b ./ se
    p  = 2 .* (1 .- cdf.(Normal(0, 1), abs.(z)))
    ci = confint(r; level=level)
    α  = 1 - level
    cl = string(round(100*(α/2);    digits=1), "%")
    cu = string(round(100*(1-α/2);  digits=1), "%")
    return CoefTable(
        hcat(b, se, z, p, ci[:, 1], ci[:, 2]),
        ["Coef.", "Std. Error", "z", "Pr(>|z|)", cl, cu],
        copy(r.coefnames_vec),
        4,  # pvalcol
        3   # teststatcol
    )
end

function Base.show(io::IO, ::MIME"text/plain", r::IVPPMLResult)
    println(io, "IVPPMLResult: ", r.converged ? "converged" : "did NOT converge",
            "  (iterations=", r.iterations, ", N=", r.N, ", N_full=", r.N_full, ")")
    println(io, "Response: ", r.depvar)
    println(io, "Deviance=", round(r.deviance; sigdigits=6),
            "   ll=", round(r.ll; sigdigits=6),
            "   ll_0=", round(r.ll_0; sigdigits=6),
            "   pseudo R²=", round(1 - r.ll / r.ll_0; digits=4))
    if r.num_sep_mu > 0
        println(io, "(", r.num_sep_mu, " obs dropped by mu-separation)")
    end
    println(io, coeftable(r))
end

Base.show(io::IO, r::IVPPMLResult) = show(io, MIME"text/plain"(), r)

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
# _select_not_collinear — rank detection on FE-partialled columns, with
# an absorption test against pre-partial (raw) columns. Mirrors ppmlhdfe's
# select_not_collinear + reghdfe_rmcoll combo and Stata ivppmlhdfe.ado's
# two-stage collinearity pipeline (v0.9.4, line 1405).
#
# Two-part test per column j:
#   (a) Absorption:  weighted norm of the partialled column must exceed
#       rtol * max(raw weighted norms). Catches columns fully absorbed by
#       the FE structure (e.g. user-expanded factor dummies matching an
#       absorbed FE — T3.1 in 02_stress_julia v0.9.2).
#   (b) Within-block rank: on the survivors, pivoted QR with a relative
#       threshold rtol * max(|R_ii|). Catches within-block collinearity
#       (e.g. z_a == z_b).
#
# Input:
#   A_partial — N×K matrix of FE-partialled columns
#   A_raw     — N×K matrix of pre-partial columns (same row filtering)
#   w         — length-N non-negative weight vector
#   rtol      — relative tolerance (default 1e-9)
#
# Output: Vector{Int} of surviving column indices (1..K), sorted.
# ============================================================================
function _select_not_collinear(A_partial::AbstractMatrix{Float64},
                               A_raw::AbstractMatrix{Float64},
                               w::AbstractVector{Float64};
                               rtol::Float64 = 1e-9)
    N, K = size(A_partial)
    K == 0 && return Int[]
    @assert size(A_raw, 2) == K "raw / partial column count mismatch"

    sw = sqrt.(max.(w, 0.0))

    # Raw-norm scale reference (absolute floor anchor)
    raw_norms = Float64[norm(sw .* view(A_raw, :, j)) for j in 1:K]
    max_raw = isempty(raw_norms) ? 0.0 : maximum(raw_norms)
    if max_raw == 0.0
        return Int[]
    end
    abs_thresh = rtol * max_raw

    # Part (a): drop columns whose partialled norm is below abs_thresh.
    partial_norms = Float64[norm(sw .* view(A_partial, :, j)) for j in 1:K]
    nonzero_idx = Int[j for j in 1:K if partial_norms[j] > abs_thresh]

    if isempty(nonzero_idx)
        return Int[]
    end
    if length(nonzero_idx) == 1
        return nonzero_idx
    end

    # Part (b): pivoted QR on non-absorbed survivors for within-block rank
    WA_sub = Matrix{Float64}(undef, N, length(nonzero_idx))
    @inbounds for (i, j) in enumerate(nonzero_idx)
        @views WA_sub[:, i] .= sw .* A_partial[:, j]
    end
    F = qr(WA_sub, ColumnNorm())
    diagR = abs.(diag(F.R))
    max_diag = isempty(diagR) ? 0.0 : maximum(diagR)
    thresh = rtol * max_diag
    kept_pivoted = findall(d -> d > thresh, diagR)
    kept_in_sub = Int[F.p[i] for i in kept_pivoted]
    kept_indices = sort!(Int[nonzero_idx[i] for i in kept_in_sub])
    return kept_indices
end

# ============================================================================
# Robust sandwich VCE  (ppmlhdfe convention: (N/(N-1)) multiplier)
# ============================================================================

function _ivppml_robust_vce(Xhat::Matrix{Float64},
                            X_dm::Matrix{Float64},
                            w::Vector{Float64},
                            resid::Vector{Float64},
                            N::Int)
    bread  = inv(Symmetric(Xhat' * (w .* X_dm)))
    scores = Xhat .* (w .* resid)
    meat   = scores' * scores
    V = (N / (N - 1)) .* bread * meat * bread
    return V
end

# ============================================================================
# Cluster meat: sum of outer products of score sums within clusters
# Returns (meat, G) where G is number of clusters
# ============================================================================

function _ivppml_clustmeat(scores::Matrix{Float64},
                           clust_id::AbstractVector)
    K = size(scores, 2)
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

    return meat, G
end

# ============================================================================
# Interaction group IDs from two group vectors
# ============================================================================

function _interact_id(id1::AbstractVector, id2::AbstractVector)
    N = length(id1)
    perm = sortperm(collect(zip(id1, id2)))
    s1 = id1[perm]
    s2 = id2[perm]
    result = similar(id1, Int)
    result[perm[1]] = 1
    cnt = 1
    for i in 2:N
        if s1[i] != s1[i-1] || s2[i] != s2[i-1]
            cnt += 1
        end
        result[perm[i]] = cnt
    end
    return result
end

# ============================================================================
# Single-cluster VCE  (Arellano G/(G-1) correction)
# ============================================================================

function _ivppml_cluster_vce_single(Xhat::Matrix{Float64},
                                    X_dm::Matrix{Float64},
                                    w::Vector{Float64},
                                    resid::Vector{Float64},
                                    clust_id::AbstractVector)
    K = size(X_dm, 2)
    bread  = inv(Symmetric(Xhat' * (w .* X_dm)))
    scores = Xhat .* (w .* resid)

    meat, G = _ivppml_clustmeat(scores, clust_id)
    if G <= 1
        # Fall back to robust VCE if only 1 cluster
        V = _ivppml_robust_vce(Xhat, X_dm, w, resid, size(X_dm, 1))
    else
        V = (G / (G - 1)) .* bread * meat * bread
    end
    return V, [G]
end

# ============================================================================
# Multi-way cluster VCE (CGM inclusion-exclusion formula)
# ============================================================================

function _ivppml_cluster_vce_multiway(Xhat::Matrix{Float64},
                                      X_dm::Matrix{Float64},
                                      w::Vector{Float64},
                                      resid::Vector{Float64},
                                      clust_ids::Matrix{Float64},
                                      n_clust::Int)
    K = size(X_dm, 2)
    bread  = inv(Symmetric(Xhat' * (w .* X_dm)))
    scores = Xhat .* (w .* resid)

    V_slope = zeros(K, K)
    G_counts = zeros(Int, n_clust)

    # Iterate over all non-empty subsets of {1..n_clust}
    for mask in 1:(2^n_clust - 1)
        subset_size = count_ones(mask)

        # Build combined cluster ID for this subset
        cid = nothing
        for j in 1:n_clust
            if (mask >> (j - 1)) & 1 == 1
                col_j = clust_ids[:, j]
                if cid === nothing
                    cid = col_j
                else
                    cid = Float64.(_interact_id(cid, col_j))
                end
            end
        end

        meat_sub, G_sub = _ivppml_clustmeat(scores, cid)

        # CGM sign: + for odd subset size, - for even
        sign = iseven(subset_size) ? -1 : 1

        if G_sub > 1
            V_slope .+= sign * (G_sub / (G_sub - 1)) .* bread * meat_sub * bread
        end

        # Store cluster counts for single-variable subsets
        if subset_size == 1
            for j in 1:n_clust
                if mask == (1 << (j - 1))
                    G_counts[j] = G_sub
                end
            end
        end
    end

    # Symmetrize
    V_slope .= 0.5 .* (V_slope .+ V_slope')

    # PSD fix for multi-way clustering (matching reghdfe)
    if n_clust > 1
        ev = eigvals(Symmetric(V_slope))
        min_ev = minimum(ev)
        if min_ev < 0
            V_slope .-= min_ev .* Matrix{Float64}(I, K, K)
        end
    end

    return V_slope, G_counts
end

# ============================================================================
# Main estimation: ivppml_reg
#
# Single-call entry point mirroring nlreg() from GLFixedEffectModels.jl.
# Accepts a DataFrame and column symbols.
# ============================================================================

# Convenience: kwargs-only entry point so users can write
#   ivppml_reg(df=df, depvar=:y, exog=[:x], endog=[:xi], instruments=[:z], fe=[:id])
# matching the FixedEffectModels.jl style.
function ivppml_reg(; df,
                    depvar::Symbol,
                    exog::AbstractVector = Symbol[],
                    endog::AbstractVector = Symbol[],
                    instruments::AbstractVector = Symbol[],
                    fe::AbstractVector = Symbol[],
                    kwargs...)
    return ivppml_reg(df, depvar, exog, endog, instruments, fe; kwargs...)
end

function ivppml_reg(df::Any,          # AbstractDataFrame
                    depvar::Symbol,
                    exog::AbstractVector,
                    endog::AbstractVector,
                    instruments::AbstractVector,
                    fe_syms::AbstractVector;
                    cluster::Union{Vector{Symbol}, Symbol, Nothing} = nothing,
                    weights::Union{Symbol, Nothing} = nothing,
                    offset::Union{Symbol, Nothing} = nothing,
                    tol::Float64  = 1e-8,
                    itol::Float64 = -1.0,
                    maxiter::Int  = 1000,
                    verbose::Int  = 0,
                    method::Symbol = :cpu,
                    standardize::Bool = false,
                    guess::Symbol = :default,
                    mu_separation::Bool = true)

    # ---- Coerce loose Vector{Any} (e.g. literal `[]`) to Vector{Symbol} ----
    # so users can write `ivppml_reg(df, :y, [:x], [], [], [:fe])` without
    # MethodError on the empty literal.
    exog        = Symbol[Symbol(s) for s in exog]
    endog       = Symbol[Symbol(s) for s in endog]
    instruments = Symbol[Symbol(s) for s in instruments]
    fe_syms     = Symbol[Symbol(s) for s in fe_syms]

    # ---- Validate keyword arguments ----
    guess in (:default, :simple, :mean) ||
        error("ivppml_reg: guess must be :default, :simple, or :mean (got :$(guess))")

    # ---- Helper: load a numeric column with informative errors ----
    function _load_col(df, sym::Symbol, role::String)
        sym in propertynames(df) ||
            error("ivppml_reg: $(role) column :$(sym) not found in DataFrame")
        col = df[!, sym]
        any(ismissing, col) &&
            error("ivppml_reg: $(role) column :$(sym) contains missing values; drop or impute before estimation")
        try
            return Float64.(col)
        catch e
            error("ivppml_reg: $(role) column :$(sym) has eltype $(eltype(col)) which cannot be converted to Float64")
        end
    end

    # Type guard: df must look like an AbstractDataFrame (or compatible).
    # df::Any in the signature invites duck typing; reject obvious non-table
    # inputs early so users see a clear message instead of a MethodError from
    # `size(::NamedTuple, ::Int)`.
    if !(applicable(size, df, 1) && applicable(propertynames, df))
        error("ivppml_reg: df must be an AbstractDataFrame-like table (got $(typeof(df)))")
    end

    # ---- Extract data from DataFrame ----
    N_full = size(df, 1)
    y = _load_col(df, depvar, "depvar")
    N = length(y)

    # Guard against pathological inputs early.
    any(!isfinite, y) && error("ivppml_reg: depvar :$(depvar) contains NaN or Inf")
    any(y .< 0)      && error("ivppml_reg: depvar :$(depvar) must be non-negative")
    sum(y) == 0      && error("ivppml_reg: depvar :$(depvar) is identically zero; model is degenerate")

    K_exog  = length(exog)
    K_endo  = length(endog)
    K       = K_exog + K_endo
    n_inst  = length(instruments)
    L       = K_exog + n_inst         # total instruments

    # Order condition: at least one instrument per endogenous regressor.
    # Without this guard the first-stage 2SLS solve `ZwZ \ ZwX` raises
    # SingularException(2) at iter 1, which the second-stage try/catch
    # does not intercept.
    n_inst < K_endo &&
        error("ivppml_reg: equation not identified — need at least $(K_endo) " *
              "instrument(s), got $(n_inst)")

    # Build X = [exog endog] using _load_col so missing/Inf/type errors are
    # informative.
    X = Matrix{Float64}(undef, N, K)
    for (j, s) in enumerate(exog)
        col = _load_col(df, s, "exogenous regressor")
        any(!isfinite, col) && error("ivppml_reg: exogenous regressor :$(s) contains NaN or Inf")
        X[:, j] .= col
    end
    for (j, s) in enumerate(endog)
        col = _load_col(df, s, "endogenous regressor")
        any(!isfinite, col) && error("ivppml_reg: endogenous regressor :$(s) contains NaN or Inf")
        X[:, K_exog + j] .= col
    end

    # Build excluded instruments
    Z_excl = Matrix{Float64}(undef, N, n_inst)
    for (j, s) in enumerate(instruments)
        col = _load_col(df, s, "instrument")
        any(!isfinite, col) && error("ivppml_reg: instrument :$(s) contains NaN or Inf")
        Z_excl[:, j] .= col
    end

    # User weights
    w_user_full = weights === nothing ? ones(Float64, N) : _load_col(df, weights, "weight")

    # Offset
    has_offset = offset !== nothing
    offset_full = if has_offset
        col = _load_col(df, offset, "offset")
        any(!isfinite, col) && error("ivppml_reg: offset :$(offset) contains NaN or Inf")
        col
    else
        zeros(Float64, N)
    end

    # Cluster IDs — normalize to Vector{Symbol}
    local clust_syms::Vector{Symbol}
    if cluster === nothing
        clust_syms = Symbol[]
    elseif cluster isa Symbol
        clust_syms = [cluster]
    else
        clust_syms = cluster
    end
    n_clust = length(clust_syms)
    has_cluster = n_clust > 0

    # Load cluster ID columns (as a Matrix with Float64 IDs, unified
    # representation so we can row-subset below).
    clust_id_mat = nothing
    if has_cluster
        clust_id_mat = Matrix{Float64}(undef, N, n_clust)
        for (j, s) in enumerate(clust_syms)
            s in propertynames(df) ||
                error("ivppml_reg: cluster column :$(s) not found in DataFrame")
            col = df[!, s]
            any(ismissing, col) &&
                error("ivppml_reg: cluster column :$(s) contains missing values")
            try
                clust_id_mat[:, j] .= Float64.(col)
            catch
                error("ivppml_reg: cluster column :$(s) has eltype $(eltype(col)) which cannot be converted to Float64 (use a numeric ID column or `groupindices(groupby(df, :$(s)))`)")
            end
        end
    end

    # Build FE raw-level vectors (need these to re-subset under weight=0 drop).
    # Reject missings here too — `FixedEffect` would otherwise propagate them
    # as NaN through demeaning and silently corrupt coefficients.
    #
    # When the user passes NO fixed effects, inject a single constant FE
    # (all ones). Without this, `solve_residuals!` is a no-op and the 2SLS
    # runs on un-centred data without an intercept column, silently producing
    # biased coefficients and pinning `b_cons = 0`.
    if isempty(fe_syms)
        fe_raw = Any[ones(Int, N_full)]
    else
        for s in fe_syms
            s in propertynames(df) ||
                error("ivppml_reg: FE column :$(s) not found in DataFrame")
            any(ismissing, df[!, s]) &&
                error("ivppml_reg: FE column :$(s) contains missing values; drop or impute before estimation")
        end
        fe_raw = Any[Vector(df[!, s]) for s in fe_syms]
    end

    # Drop obs with zero or negative weight BEFORE anything else — matches
    # Stata's behaviour with [pw=...]. Otherwise NaN propagates into the
    # 2SLS solve via getrf! on a rank-deficient design.
    # Track which rows survived so r.esample correctly reflects the
    # ESTIMATION sample (not the original df row count).
    keep_full = trues(N_full)
    if weights !== nothing
        any(!isfinite, w_user_full) && error("ivppml_reg: weights contain NaN/Inf")
        any(w_user_full .< 0)        && error("ivppml_reg: weights must be non-negative")
        keep = w_user_full .> 0
        if !all(keep)
            y           = y[keep]
            X           = X[keep, :]
            Z_excl      = Z_excl[keep, :]
            w_user_full = w_user_full[keep]
            offset_full = offset_full[keep]
            if has_cluster
                clust_id_mat = clust_id_mat[keep, :]
            end
            fe_raw = Any[v[keep] for v in fe_raw]
            keep_full = keep
            N = length(y)
        end
    end

    w_user     = w_user_full
    offset_vec_orig = offset_full
    # Center the offset to weighted mean 0 for the IRLS loop.  Without this a
    # large absolute-scale offset (e.g. log(trade_$) ≈ 25) leaves the initial
    # mu = exp(offset) ≈ 1e11 and explodes the bread matrix.  After the loop
    # we restore `offset_vec` to the user's original values for b_cons / ll_0
    # / d_vals — eta / mu / irls_w are already on the user's scale because
    # the IRLS absorbed the centering shift into the fitted b_cons.  Mirrors
    # native ivppmlhdfe and ppmlhdfe.mata:461.
    offset_mean = has_offset ? dot(w_user, offset_full) / sum(w_user) : 0.0
    offset_vec  = has_offset ? (offset_full .- offset_mean) : offset_full

    # Build FE objects from the (possibly row-filtered) raw vectors.
    fes = [FixedEffect(v) for v in fe_raw]

    # Expose cluster matrix back to downstream code under its old name.
    clust_id_cols = nothing
    if has_cluster
        clust_id_cols = n_clust == 1 ? clust_id_mat[:, 1] : clust_id_mat
    end

    # Coefficient names
    cnames = vcat(String.(vcat(exog, endog)), "_cons")

    # ---- Pre-IRLS sanity check on X / Z columns ----
    # Reject zero-variance regressors and identical columns BEFORE IRLS, so
    # the user gets a clear error rather than 1e12 garbage coefs / a singular
    # 2SLS solve at iter 1.
    if K > 0
        sw_check = sum(w_user)
        @inbounds for k in 1:K
            mx = dot(w_user, view(X, :, k)) / sw_check
            vx = dot(w_user, (view(X, :, k) .- mx) .^ 2) / sw_check
            if vx <= 0 || !isfinite(vx)
                varname = k <= K_exog ? exog[k] : endog[k - K_exog]
                error("ivppml_reg: regressor :$(varname) has zero variance " *
                      "(constant column); drop it before estimation")
            end
        end
    end
    if n_inst > 0
        sw_check = sum(w_user)
        @inbounds for k in 1:n_inst
            mz = dot(w_user, view(Z_excl, :, k)) / sw_check
            vz = dot(w_user, (view(Z_excl, :, k) .- mz) .^ 2) / sw_check
            if vz <= 0 || !isfinite(vz)
                error("ivppml_reg: instrument :$(instruments[k]) has zero " *
                      "variance (constant column); drop it before estimation")
            end
        end
    end

    # ---- Standardize X, Z (optional) ----
    # Divides X, Z by column stdev for numerical stability.
    # Coefficients and VCE are back-transformed after convergence.
    stdev_x = ones(Float64, K)
    stdev_z = ones(Float64, n_inst)
    if standardize && K > 0
        sw_std = sum(w_user)
        for k in 1:K
            mx = dot(w_user, view(X, :, k)) / sw_std
            v  = dot(w_user, (view(X, :, k) .- mx) .^ 2) / sw_std
            sd = sqrt(max(v, 0.0))
            if sd > 0
                stdev_x[k] = sd
                X[:, k] ./= sd
            end
        end
        for k in 1:n_inst
            mz = dot(w_user, view(Z_excl, :, k)) / sw_std
            v  = dot(w_user, (view(Z_excl, :, k) .- mz) .^ 2) / sw_std
            sd = sqrt(max(v, 0.0))
            if sd > 0
                stdev_z[k] = sd
                Z_excl[:, k] ./= sd
            end
        end
    end

    # Per-block collinearity drop counts (propagated to IVPPMLResult for bridge).
    n_drop_exog  = 0
    n_drop_endog = 0
    n_drop_inst  = 0

    # ---- Two-stage collinearity removal (mirrors Stata ivppmlhdfe v0.9.4) ----
    # Detect linearly-dependent columns within exog, endog, and instrument
    # blocks AFTER projecting them through the FE structure. Three independent
    # passes (one per IV block) because IV-PPML needs each matrix full-rank
    # for the 2SLS solve to be well-posed. Without this, factor-variable
    # expansions like manually-constructed `i.x + absorb(x)` or ibn-collinear
    # instruments silently produced rank-deficient designs that LU-pivoted
    # `\` silently returned garbage for (see T3.1 "manual xf dummies" in the
    # v0.9.2 stress log — coefs of order 1e13 with converged=false).
    if (K + n_inst) > 0
        target_inner_tol_coll = max(1e-12, 0.1 * tol)
        feM_coll = AbstractFixedEffectSolver{Float64}(fes, Weights(w_user), Val{method})
        # Keep a RAW (un-partialled) copy for the absorption test in
        # _select_not_collinear — we need both matrices to detect columns
        # whose signal was fully absorbed by the FE structure.
        data_raw = Matrix{Float64}(undef, N, K + n_inst)
        if K > 0
            data_raw[:, 1:K] .= X
        end
        if n_inst > 0
            data_raw[:, K+1:end] .= Z_excl
        end
        data_coll = copy(data_raw)
        solve_residuals!(data_coll, feM_coll;
                         maxiter=10000, tol=target_inner_tol_coll, progress_bar=false)

        exog_keep  = K_exog > 0 ?
            _select_not_collinear(view(data_coll, :, 1:K_exog),
                                  view(data_raw,  :, 1:K_exog),
                                  w_user) : Int[]
        endog_keep = K_endo > 0 ?
            _select_not_collinear(view(data_coll, :, (K_exog+1):K),
                                  view(data_raw,  :, (K_exog+1):K),
                                  w_user) : Int[]
        inst_keep  = n_inst > 0 ?
            _select_not_collinear(view(data_coll, :, (K+1):(K+n_inst)),
                                  view(data_raw,  :, (K+1):(K+n_inst)),
                                  w_user) : Int[]

        n_drop_exog  = K_exog - length(exog_keep)
        n_drop_endog = K_endo - length(endog_keep)
        n_drop_inst  = n_inst - length(inst_keep)

        if n_drop_exog + n_drop_endog + n_drop_inst > 0
            if verbose > -1
                @printf("(two-stage collinearity: dropped %d exog, %d endog, %d inst)\n",
                        n_drop_exog, n_drop_endog, n_drop_inst)
            end

            # Re-slice X / stdev_x — exog come first, then endog
            e_idx = collect(exog_keep)
            d_idx = Int[i + K_exog for i in endog_keep]
            x_idx = vcat(e_idx, d_idx)
            if !isempty(x_idx)
                X       = X[:, x_idx]
                stdev_x = stdev_x[x_idx]
            else
                X       = Matrix{Float64}(undef, N, 0)
                stdev_x = Float64[]
            end

            # Re-slice Z_excl / stdev_z
            if n_inst > 0
                Z_excl  = Z_excl[:, inst_keep]
                stdev_z = stdev_z[inst_keep]
            end

            # Rebuild Symbol vectors and counts
            exog        = Symbol[exog[i]        for i in exog_keep]
            endog       = Symbol[endog[i]       for i in endog_keep]
            instruments = Symbol[instruments[i] for i in inst_keep]
            K_exog = length(exog_keep)
            K_endo = length(endog_keep)
            K      = K_exog + K_endo
            n_inst = length(inst_keep)
            L      = K_exog + n_inst
            cnames = vcat(String.(vcat(exog, endog)), "_cons")

            # Order condition (post-drop): re-check identification
            n_inst < K_endo &&
                error("ivppml_reg: collinearity removal left equation under-identified — " *
                      "$(K_endo) endog regressor(s) need at least $(K_endo) instrument(s), " *
                      "got $(n_inst)")
        end

        feM_coll = nothing
    end

    # ---- Initialise mu ----
    sw = sum(w_user)
    mean_y = dot(w_user, y) / sw
    if guess === :mean
        mu = fill(mean_y, N)
    else
        mu = @. 0.5 * (y + mean_y)
    end
    mu  = max.(mu, max.(0.05 .* y, 1e-3))
    eta = log.(mu) .+ offset_vec

    # Mu-separation mask (accumulates across iterations)
    sep_mask = falses(N)
    num_sep_mu = 0

    # ---- Build FE solver ONCE ----
    irls_w = w_user .* mu
    irls_w = max.(irls_w, 1e-20)
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
    target_inner_tol = itol > 0 ? itol : max(1e-12, 0.1 * tol)
    inner_tol = max(start_inner_tol, tol)
    alt_tol   = start_inner_tol

    # Step-halving
    max_step_halving   = 2
    step_memory        = 0.9
    num_step_halving   = 0

    # Pre-allocate
    b     = zeros(K)
    b_old = zeros(K)
    resid = zeros(N)
    Xhat  = zeros(N, K)
    old_eta = similar(eta)
    beta_change = Inf

    ncols = 1 + K + n_inst
    data  = Matrix{Float64}(undef, N, ncols)
    X_dm  = Matrix{Float64}(undef, N, K)
    Z_dm  = Matrix{Float64}(undef, N, L)

    # First-stage F diagnostic emitted once, after iter 1.  K_endo > 1 is not
    # implemented here (that needs Cragg-Donald or Kleibergen-Paap).
    fstat_warned = false

    if verbose > -1
        @printf("IRLS-IV iterations (N = %d, K = %d, L = %d)\n", N, K, L)
        println("-"^60)
    end

    for it in 1:maxiter
        iter = it

        # (a0) Early divergence guard — matches Stata ivppmlhdfe v0.9.3 rc=9003.
        # Catches Class C pathology where eta diverges before solve_residuals!
        # or 2SLS (z = eta - 1 + y/mu becomes non-finite), and runaway |β|.
        if any(!isfinite, mu) || any(!isfinite, eta)
            error("ivppml_reg: mu/eta has Inf/NaN at iter $(iter) (rc=9003 runaway " *
                  "divergence); try mu_separation=true, standardize=true, or simpler FE")
        end
        if iter > 10 && K > 0 &&
           (maximum(abs, b) > 1e6 || maximum(abs, eta) > 30)
            msg = @sprintf("ivppml_reg: IRLS diverged at iter %d (rc=9003: max|b|=%.2e, max|eta|=%.2e); aborting",
                           iter, maximum(abs, b), maximum(abs, eta))
            error(msg)
        end

        # (a) Working dependent variable: z = eta - offset - 1 + y / mu
        z = eta .- offset_vec .- 1.0 .+ y ./ mu
        # For y=0 obs under mu-separation: pin z = eta - offset - 1 exactly
        # (avoid 0/epsilon noise, matching ppmlhdfe zero_sample handling)
        if mu_separation
            @inbounds for i in 1:N
                if y[i] == 0.0
                    z[i] = eta[i] - offset_vec[i] - 1.0
                end
            end
        end

        # (b) IRLS weights: w = w_user * mu (clamp to keep solver alive on overshoot)
        @. irls_w = w_user * mu
        @. irls_w = clamp(irls_w, 1e-20, 1e20)

        # Runaway-divergence guard — matches Stata ivppmlhdfe v0.9.3 rc=9003
        # (and ppmlhdfe.mata line 620). If mu/irls_w goes Inf/NaN (typically
        # because eta blew up under repeated step-halving failures), hard-error
        # rather than leaking junk results downstream.
        if any(!isfinite, irls_w) || any(!isfinite, mu)
            error("ivppml_reg: mu has Inf/NaN at iter $(iter) (rc=9003 runaway " *
                  "divergence); try mu_separation=true, standardize=true, or simpler FE")
        end

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

        # (e) Weighted 2SLS, with collinearity guard.
        # The 2SLS \ solve can throw SingularException on perfectly collinear
        # X. Catch it, surface a clear error, and stop (rather than leaking a
        # raw stacktrace into the caller).
        if K > 0; copyto!(b_old, b); end
        try
            _ivppml_2sls!(b, resid, Xhat, z_dm, X_dm, Z_dm, irls_w)
        catch e
            if e isa LinearAlgebra.SingularException
                error("ivppml_reg: 2SLS solve failed at iter $(iter) — the " *
                      "design matrix appears rank-deficient. Check for " *
                      "perfectly collinear regressors / instruments, or set " *
                      "standardize=true.")
            else
                rethrow(e)
            end
        end
        if K > 0; beta_change = maximum(abs.(b .- b_old)); end

        # (e2) First-stage weak-IV F warning (one-shot at iter 1, K_endo==1).
        # Partial F-stat for joint significance of excluded instruments after
        # partialling out the exogenous regressors. Staiger-Stock threshold 10.
        if !fstat_warned && iter == 1 && K_endo == 1 && n_inst > 0
            y_fs = view(X_dm, :, K)   # single endog column sits last in X_dm
            # Full model: OLS y_fs on Z_dm (exog + excl instruments)
            ZwZ_fs = Z_dm' * (irls_w .* Z_dm)
            Zwy_fs = Z_dm' * (irls_w .* y_fs)
            beta_full = ZwZ_fs \ Zwy_fs
            y_hat_full = Z_dm * beta_full
            ssr_full = sum(irls_w .* (y_fs .- y_hat_full).^2)
            # Restricted: OLS y_fs on exog only (or 0-mean if K_exog==0)
            ssr_r = 0.0
            if K_exog > 0
                X_exog_dm = view(X_dm, :, 1:K_exog)
                XwX = X_exog_dm' * (irls_w .* X_exog_dm)
                Xwy = X_exog_dm' * (irls_w .* y_fs)
                beta_r = XwX \ Xwy
                y_hat_r = X_exog_dm * beta_r
                ssr_r = sum(irls_w .* (y_fs .- y_hat_r).^2)
            else
                ssr_r = sum(irls_w .* (y_fs .^ 2))
            end
            df_num = n_inst
            df_den = max(N - L, 1)
            if ssr_full > 0 && df_num > 0
                F_fs = (ssr_r - ssr_full) / df_num / (ssr_full / df_den)
                if F_fs < 10
                    @warn "first-stage F < 10 (weak instruments — inference may be unreliable)" F_fs
                end
            end
            fstat_warned = true
        end

        # (f) Update eta = z - resid + offset
        copyto!(old_eta, eta)
        @. eta = z - resid + offset_vec

        # (f2) Mu-separation check (matching ppmlhdfe).
        # Detects y=0 obs whose eta is diverging to -infinity.
        # Once flagged, stays flagged for the rest of estimation.
        if mu_separation && iter > 1
            log_septol = log(1e-6)
            # Track the SMALLEST eta among y>0 obs (the "still-active" subspace).
            # The original v0.9.2 code initialized min_eta_pos = -Inf and updated
            # only when eta[i] < min_eta_pos, which is impossible — so the
            # adaptive branch never fired. Initialize with +Inf instead.
            min_eta_pos = Inf
            @inbounds for i in 1:N
                if y[i] > 0.0 && eta[i] < min_eta_pos
                    min_eta_pos = eta[i]
                end
            end
            adjusted_log_septol = isfinite(min_eta_pos) ?
                log_septol + min(min_eta_pos + 5.0, 0.0) : log_septol
            n_new_sep = 0
            @inbounds for i in 1:N
                if !sep_mask[i] && y[i] == 0.0 && eta[i] <= adjusted_log_septol
                    sep_mask[i] = true
                    n_new_sep += 1
                end
            end
            if n_new_sep > 0
                num_sep_mu += n_new_sep
                if verbose > -1
                    @printf("(mu-separation: %d obs detected at iter %d)\n", n_new_sep, iter)
                end
            end
        end

        # (g) Update mu
        @. mu = exp(eta)
        # Set mu=0 for separated obs (matching ppmlhdfe), then floor
        if mu_separation && num_sep_mu > 0
            @inbounds for i in 1:N
                if sep_mask[i]
                    mu[i] = 0.0
                end
            end
        end
        @. mu = max(mu, eps(100.0))

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

        # Numerical safeguards (matching ppmlhdfe edittozerotol)
        if 2.0 * deviance / N < eps(1.0)
            deviance = 0.0
        end

        # (i) Convergence + step-halving + adaptive tolerance
        is_step_halving = false
        if iter > 1
            delta_dev = old_deviance - deviance
            # Clip delta_dev if deviance is small relative to delta
            if deviance < 0.1 * delta_dev
                delta_dev = deviance
            end
            denom     = max(min(deviance, old_deviance), 0.1)
            eps_val   = abs(delta_dev) / denom

            if eps_val < tol
                if inner_tol <= 1.1 * target_inner_tol || length(fes) <= 1
                    ok += 1
                    ok >= 1 && (converged = true)
                end
            elseif delta_dev < 0 && num_step_halving < max_step_halving
                @. eta = step_memory * old_eta + (1 - step_memory) * eta
                # Clip eta at -10 on repeated step-halving (matching ppmlhdfe)
                if num_step_halving > 0
                    @. eta = max(eta, -10.0)
                end
                @. mu  = exp(eta)
                @. mu  = max(mu, eps(100.0))
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
            K > 0 && isfinite(beta_change) && @printf("  db = %-9.4e", beta_change)
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

    if !converged
        # Matches Stata ivppmlhdfe v0.9.3 rc=430 (loud-failure on non-convergence).
        msg = @sprintf("ivppml_reg: failed to converge in %d iterations (rc=430 non-convergence, eps=%.4e); try increasing maxiter, or use standardize=true / mu_separation=true",
                       maxiter, eps_val)
        error(msg)
    elseif verbose > -1
        @printf("Converged in %d iterations (tol = %9.4e)\n", iter, tol)
    end

    # ==================================================================
    # Final beta with constant, and VCE
    # ==================================================================

    @. irls_w = max(w_user * mu, 1e-20)
    sw_f = sum(irls_w)

    # Restore offset to the user's original (uncentered) values.  eta and mu
    # are already on the user's scale (the IRLS absorbed the centring shift
    # into the fitted intercept), so only `offset_vec` needs restoring for
    # the post-processing formulas below (b_cons / ll_0 / d_vals).
    if has_offset && offset_mean != 0.0
        offset_vec = offset_vec_orig
    end

    mean_eta_no_offset = dot(irls_w, eta .- offset_vec) / sw_f
    mean_X   = vec(sum(irls_w .* X; dims=1)) ./ sw_f
    b_cons   = mean_eta_no_offset - dot(mean_X, b)

    b_full  = vcat(b, b_cons)
    K_total = K + 1

    # Recompute Xhat for VCE
    ZwZ_f = Z_dm' * (irls_w .* Z_dm)
    Pi_f  = ZwZ_f \ (Z_dm' * (irls_w .* X_dm))
    Xhat  .= Z_dm * Pi_f

    # Rank check (collinearity): rank of bread matrix Xhat'WX
    actual_rank = K
    if K > 0
        bread_mat = Xhat' * (irls_w .* X_dm)
        actual_rank = rank(bread_mat; rtol = 1e-12)
    end

    # VCE
    G_counts = Int[]
    if has_cluster
        if n_clust == 1
            V_slope, G_counts = _ivppml_cluster_vce_single(Xhat, X_dm, irls_w, resid, clust_id_cols)
        else
            V_slope, G_counts = _ivppml_cluster_vce_multiway(Xhat, X_dm, irls_w, resid, clust_id_cols, n_clust)
        end
    else
        V_slope = _ivppml_robust_vce(Xhat, X_dm, irls_w, resid, N)
    end

    # ---- Back-transform standardized coefficients and VCE ----
    # b_cons is invariant to standardization (mean(X_std)*b_std = mean(X)*b_raw),
    # so we recompute b_full after rescaling b.
    if standardize && K > 0
        b ./= stdev_x
        V_slope ./= (stdev_x * stdev_x')
        # Restore X to original scale so d_vals computation is correct below
        for k in 1:K
            X[:, k] .*= stdev_x[k]
        end
        # Restore Z_excl as well (kept for completeness; not used downstream)
        for k in 1:n_inst
            Z_excl[:, k] .*= stdev_z[k]
        end
        # Recompute b_cons + b_full on the un-standardized scale
        mean_X_orig = vec(sum(irls_w .* X; dims=1)) ./ sw_f
        b_cons = mean_eta_no_offset - dot(mean_X_orig, b)
        b_full = vcat(b, b_cons)
    end

    # Expand to include constant (VCE for _cons is zero: not estimable with absorbed FE)
    # Following ivreg2 partial() convention
    V = zeros(K_total, K_total)
    if K > 0
        V[1:K, 1:K] = V_slope
    end

    # ---- Log pseudo-likelihood ----
    # ll = sum(w * (y * eta - mu - lngamma(y + 1)))
    ll = 0.0
    @inbounds for i in 1:N
        ll += w_user[i] * (y[i] * eta[i] - mu[i] - loggamma(y[i] + 1.0))
    end

    # Null model: intercept-only, offset-aware
    ll_0 = 0.0
    if has_offset && sum(abs.(offset_vec)) > 0.0
        # With offset: mu_0i = exp(c + o_i), c = log(sum(w*y) / sum(w*exp(o)))
        ll_0_c = log(dot(w_user, y) / dot(w_user, exp.(offset_vec)))
        @inbounds for i in 1:N
            mu_0i = exp(ll_0_c + offset_vec[i])
            ll_0 += w_user[i] * (y[i] * (ll_0_c + offset_vec[i]) - mu_0i - loggamma(y[i] + 1.0))
        end
    else
        # Without offset: mu_0 = mean(y, w). Guard mean_y == 0 (can occur if
        # the IRLS broke down before mean_y could be recomputed for the trimmed
        # sample) so log(0) doesn't poison ll_0 with -Inf/NaN.
        if mean_y > 0
            log_mean_y = log(mean_y)
            @inbounds for i in 1:N
                ll_0 += w_user[i] * (y[i] * log_mean_y - mean_y - loggamma(y[i] + 1.0))
            end
        else
            ll_0 = 0.0
        end
    end

    # ---- d values: FE sum = eta - offset - X*b - b_cons ----
    d_vals = eta .- offset_vec
    if K > 0
        d_vals .-= X * b
    end
    d_vals .-= b_cons
    # Center d to weighted mean zero
    d_vals .= d_vals .- dot(irls_w, d_vals) / sum(irls_w)

    # esample (all true since data is pre-filtered by Stata)
    # esample reflects which rows of the ORIGINAL df entered the estimation
    # sample. Equals all-true unless the weight=0 drop trimmed rows above.
    esample = BitVector(keep_full)

    # Absorbed-FE degrees of freedom (first-order approximation: sum of FE
    # levels minus the number of FEs, which is exact in the absence of
    # nesting/redundancy across FEs).  When fe_syms was empty, we synthesised
    # a single all-ones FE just to satisfy FixedEffects.jl, so the contribution
    # should be 0.  Subtract one per real FE to account for the intercept that
    # absorbs one level per FE.
    df_absorb = isempty(fe_syms) ? 0 :
        sum(length(unique(fe.refs)) for fe in fes) - length(fes)

    return IVPPMLResult(b_full, V, cnames, depvar, N, N_full,
                        converged, iter, deviance, ll, ll_0,
                        d_vals, eta, G_counts, esample,
                        num_sep_mu, actual_rank, df_absorb,
                        n_drop_exog, n_drop_endog, n_drop_inst)
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
Base.propertynames(d::_DFWrap) = propertynames(d.nt)

@compile_workload begin
    # Deterministic, well-conditioned dummy problem: N=40, K=1 exog, 1 endog,
    # 1 instrument, 1 FE. No random noise — converges in ~3 IRLS iters under
    # all seeds, satisfying v0.9.3 loud-failure rc=430 / rc=9003 guards.
    _N = 40
    _t = Float64.(1:_N)
    _x1 = 0.10 .* _t
    _x2 = 0.05 .* _t
    _z1 = _x2 .+ 0.01 .* _t
    _eta = 0.10 .+ 0.20 .* _x1 .+ 0.05 .* _x2
    _y = exp.(_eta)
    _dfwrap = _DFWrap((;
        y    = _y,
        x1   = _x1,
        x2   = _x2,
        z1   = _z1,
        fe1  = mod.(0:_N-1, 5) .+ 1,
        cid  = mod.(0:_N-1, 8) .+ 1,
        off1 = zeros(_N),
    ))

    # Run estimation (silent) — triggers compilation of all code paths.
    # Wrap in try/catch: if the tiny problem hits a numerical edge case the
    # precompile should still succeed; first real call will JIT-compile any
    # uncovered code path.
    try
        ivppml_reg(_dfwrap, :y, Symbol[:x1], Symbol[:x2], Symbol[:z1], Symbol[:fe1];
                   cluster = :cid, offset = :off1, tol = 1e-8, maxiter = 50, verbose = -1)
    catch
    end
end

# ============================================================================
# ppml_reg — convenience wrapper: PPML without IV
# When no endogenous regressors / instruments are provided, the IRLS-IV
# algorithm degenerates to IRLS with OLS (Z_dm = X_dm, Xhat = X_dm).
# This gives identical results to GLFixedEffectModels / ppmlhdfe.
# ============================================================================

"""
    ppml_reg(df, depvar, regressors, fe_syms; kwargs...)

Poisson PML with high-dimensional fixed effects (no IV).
Equivalent to `ivppml_reg` with empty `endog` and `instruments`.
"""
function ppml_reg(df::Any,
                  depvar::Symbol,
                  regressors::AbstractVector,
                  fe_syms::AbstractVector;
                  cluster::Union{Vector{Symbol}, Symbol, Nothing} = nothing,
                  weights::Union{Symbol, Nothing} = nothing,
                  offset::Union{Symbol, Nothing} = nothing,
                  tol::Float64  = 1e-8,
                  itol::Float64 = -1.0,
                  maxiter::Int  = 1000,
                  verbose::Int  = 0,
                  method::Symbol = :cpu,
                  standardize::Bool = false,
                  guess::Symbol = :default,
                  mu_separation::Bool = true)
    return ivppml_reg(df, depvar, regressors, Symbol[], Symbol[], fe_syms;
                      cluster=cluster, weights=weights, offset=offset,
                      tol=tol, itol=itol, maxiter=maxiter, verbose=verbose,
                      method=method, standardize=standardize, guess=guess,
                      mu_separation=mu_separation)
end

end # module
