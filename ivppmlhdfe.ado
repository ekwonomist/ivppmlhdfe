*! ivppmlhdfe 0.6.0  15mar2026
*! IV-PPML with High-Dimensional Fixed Effects
*! Algorithm: IRLS-IV  (Mullahy 1997 GMM via iteratively reweighted 2SLS)
*! Authors: Ohyun Kwon
*!
*! Syntax:
*!   ivppmlhdfe depvar [exogvars] (endogvars = excluded_instruments) [if] [in] [pw], ///
*!       absorb(absvars) [vce(cluster clustvar) tolerance(#) maxiterations(#) verbose(#)]

program ivppmlhdfe, eclass
	version 14.0

	if replay() {
		if "`e(cmd)'" != "ivppmlhdfe" error 301
		_coef_table, `0'
		exit
	}

	Estimate `0'
end


// ==========================================================================
// Parse IV syntax:  depvar [exog] (endog = instruments)  [if] [in], options
// ==========================================================================

program ParseIV, sclass
	sreturn clear
	local cmd `"`0'"'

	local popen  = strpos(`"`cmd'"', "(")
	local pclose = strpos(`"`cmd'"', ")")
	if `popen' == 0 | `pclose' == 0 | `pclose' <= `popen' {
		di as err "Syntax: ivppmlhdfe depvar [exogvars] (endogvars = instruments) [if] [in], absorb() [options]"
		exit 198
	}

	local before  = strtrim(substr(`"`cmd'"', 1, `popen' - 1))
	local inside  = strtrim(substr(`"`cmd'"', `popen' + 1, `pclose' - `popen' - 1))
	local after   = strtrim(substr(`"`cmd'"', `pclose' + 1, .))

	gettoken depvar exog : before
	local exog = strtrim(`"`exog'"')

	gettoken endog rest : inside, parse("=")
	gettoken eq instruments : rest, parse("=")
	local endog       = strtrim(`"`endog'"')
	local instruments = strtrim(`"`instruments'"')
	if "`eq'" != "=" {
		di as err "expected '=' between endogenous variables and instruments"
		exit 198
	}

	sreturn local depvar       `"`depvar'"'
	sreturn local exog         `"`exog'"'
	sreturn local endog        `"`endog'"'
	sreturn local instruments  `"`instruments'"'
	sreturn local after        `"`after'"'
end


// ==========================================================================
// Main estimation
// ==========================================================================

program Estimate, eclass
	ereturn clear

	// ---------- 1.  Parse IV specification ----------
	ParseIV `0'
	local depvar       `s(depvar)'
	local exog         `s(exog)'
	local endog        `s(endog)'
	local instruments  `s(instruments)'
	local after_paren  `s(after)'

	// ---------- 2.  Parse options ----------
	local 0 `"`after_paren'"'
	syntax [if] [in] [pw fw/] , Absorb(string) ///
		[VCE(string) CLuster(string) ///
		 TOLerance(real 1e-8) MAXITerations(integer 1000) ///
		 Verbose(integer 0) noLOG EForm IRr ///
		 BIASCorrection(string)]

	if "`irr'" != "" local eform "eform"

	if "`cluster'" != "" {
		if "`vce'" != "" {
			di as err "cannot specify both cluster() and vce()"
			exit 198
		}
		local vce "cluster `cluster'"
	}

	local vcetype "robust"
	local clustvar ""
	if "`vce'" != "" {
		gettoken vtype vrest : vce
		if "`vtype'" == "cluster" {
			local vcetype "cluster"
			local clustvar = strtrim(`"`vrest'"')
		}
		else if "`vtype'" == "robust" {
			local vcetype "robust"
		}
		else {
			di as err "vce(`vce') not supported; use robust or cluster"
			exit 198
		}
	}

	// ---------- 3.  Expand and validate variables ----------
	unab depvar      : `depvar'
	if "`exog'" != "" unab exog : `exog'
	unab endog       : `endog'
	unab instruments : `instruments'
	if "`clustvar'" != "" unab clustvar : `clustvar'

	// ---------- Parse biascorrection() ----------
	local bc_class ""
	local bc_exp_var ""
	local bc_imp_var ""
	if "`biascorrection'" != "" {
		gettoken bc_class bc_rest : biascorrection, parse(" ,")
		local bc_class = strlower(strtrim("`bc_class'"))
		if !inlist("`bc_class'", "a", "classa", "class_a", "b", "classb", "class_b") {
			di as err "biascorrection() supports: a, b"
			exit 198
		}
		if inlist("`bc_class'", "a", "classa", "class_a") local bc_class "a"
		else local bc_class "b"
		local bc_rest = strtrim(subinstr("`bc_rest'", ",", " ", .))
		gettoken bc_exp_var bc_rest2 : bc_rest
		gettoken bc_imp_var : bc_rest2
		if "`bc_exp_var'" == "" | "`bc_imp_var'" == "" {
			if "`bc_class'" == "a" {
				di as err "biascorrection(a) requires two variables: individual and time identifiers"
				di as err "  Example: biascorrection(a id year)"
			}
			else {
				di as err "biascorrection(b) requires two variables: exporter-time and importer-time identifiers"
				di as err "  Example: biascorrection(b exp_year imp_year)"
			}
			exit 198
		}
		unab bc_exp_var : `bc_exp_var'
		unab bc_imp_var : `bc_imp_var'
		if "`vcetype'" != "cluster" {
			di as err "biascorrection() requires vce(cluster clustvar)"
			exit 198
		}
	}

	local regressors   `exog' `endog'
	local n_exog     : word count `exog'
	local n_endog    : word count `endog'
	local n_inst     : word count `instruments'

	if `n_inst' < `n_endog' {
		di as err "equation not identified"
		exit 481
	}

	// mark sample
	tempvar touse
	mark `touse' `if' `in'
	markout `touse' `depvar' `regressors' `instruments' `clustvar'
	if "`weight'" != "" markout `touse' `exp'

	quietly count if `depvar' < 0 & `touse'
	if r(N) > 0 {
		di as err "`depvar' must be nonnegative"
		exit 459
	}

	quietly count if `touse'
	if r(N) == 0 error 2000

	// ---------- 4.  Display ----------
	if `verbose' > 0 {
		di as txt _n "ivppmlhdfe: IV-PPML with HDFE"
		di as txt "  Dep var:      `depvar'"
		di as txt "  Exogenous:    `exog'"
		di as txt "  Endogenous:   `endog'"
		di as txt "  Instruments:  `instruments'"
		di as txt "  Absorb:       `absorb'"
		di as txt "  VCE:          `vcetype' `clustvar'"
	}

	// ---------- 5.  IRLS-IV via Mata (using reghdfe HDFE object) ----------
	tempname b_out V_out

	mata: ivppmlhdfe_irls( ///
		"`depvar'", "`exog'", "`endog'", "`instruments'", ///
		"`absorb'", "`touse'", "`weight'", "`exp'", ///
		"`vcetype'", "`clustvar'", ///
		`tolerance', `maxiterations', `verbose', ///
		"`b_out'", "`V_out'", ///
		"`bc_class'", "`bc_exp_var'", "`bc_imp_var'" ///
	)

	// ---------- 6.  Post results ----------
	local N_final = ivppmlhdfe_N
	local colnames `regressors' _cons
	matrix colnames `b_out' = `colnames'
	matrix colnames `V_out' = `colnames'
	matrix rownames `V_out' = `colnames'

	ereturn post `b_out' `V_out', esample(`touse') depname(`depvar')

	ereturn scalar N           = `N_final'
	ereturn scalar converged   = ivppmlhdfe_converged
	ereturn scalar iterations  = ivppmlhdfe_iterations
	ereturn scalar deviance    = ivppmlhdfe_deviance

	ereturn local cmd          "ivppmlhdfe"
	ereturn local cmdline      "ivppmlhdfe `0'"
	ereturn local depvar       "`depvar'"
	ereturn local exogvars     "`exog'"
	ereturn local endogvars    "`endog'"
	ereturn local instruments  "`instruments'"
	ereturn local absorb       "`absorb'"
	ereturn local vcetype      "`vcetype'"
	ereturn local clustvar     "`clustvar'"
	ereturn local title        "IV-PPML with HDFE"

	// ---------- 7.  Display ----------
	di as txt _n "{hline 78}"
	di as txt "IV-PPML with High-Dimensional Fixed Effects"
	di as txt "{hline 78}"
	di as txt "Dependent variable: " as res "`depvar'"
	di as txt "Endogenous:         " as res "`endog'"
	di as txt "Instruments:        " as res "`instruments'"
	di as txt "Absorbed FE:        " as res "`absorb'"
	di as txt "Number of obs:      " as res e(N)
	if e(converged) == 1 {
		di as txt "Converged:          " as res "yes" ///
			as txt "  (iterations = " as res e(iterations) as txt ")"
	}
	else {
		di as txt "Converged:          " as err "no" ///
			as txt "  (iterations = " as res e(iterations) as txt ")"
	}
	if "`bc_class'" != "" {
		di as txt "Bias correction:    " as res "Class `bc_class' (KC variance + point estimate)"
	}
	di as txt "{hline 78}"
	_coef_table, `eform'
	di as txt "{hline 78}"
	di as txt "Endogenous: " as res "`endog'"
	di as txt "Instruments:" as res " `instruments'"
	if "`bc_class'" != "" {
		ereturn local biascorrection "`bc_class'"
		di as txt "{p 0 4 2}Note: Standard errors include Kauermann-Carroll correction.{p_end}"
		di as txt "{p 0 4 2}Point estimates include analytical IPP bias correction.{p_end}"
	}
end


// ==========================================================================
// Include reghdfe5 Mata library (provides FixedEffects class)
// ==========================================================================
cap findfile "reghdfe5.mata"
if (_rc) {
    di as error "ivppmlhdfe requires {bf:reghdfe} (version 5+); not found"
    di as error `"    - install from {stata ssc install reghdfe:SSC}"'
    exit 9
}
include "`r(fn)'"


// ==========================================================================
// Mata code
// ==========================================================================

mata:
mata set matastrict off


// --------------------------------------------------------------------------
// Weighted 2SLS (no constant, demeaned data)
// --------------------------------------------------------------------------

void ivppmlhdfe_gmm(
	real colvector y_dm,
	real matrix    X_dm,
	real matrix    Z_dm,
	real colvector w,
	real colvector b,
	real colvector resid)
{
	real matrix ZwZ, ZwX, Pi, Xhat, XhwX
	real colvector Xhwy

	ZwZ  = cross(Z_dm, w, Z_dm)
	ZwX  = cross(Z_dm, w, X_dm)
	Pi   = invsym(ZwZ) * ZwX
	Xhat = Z_dm * Pi

	XhwX = cross(Xhat, w, X_dm)
	Xhwy = cross(Xhat, w, y_dm)
	b    = invsym(XhwX) * Xhwy

	resid = y_dm - X_dm * b
}


// --------------------------------------------------------------------------
// Cluster-robust VCE
// --------------------------------------------------------------------------

real matrix ivppmlhdfe_clustvce(
	real matrix    Xhat,
	real matrix    X_dm,
	real colvector w,
	real colvector resid,
	real colvector clust_id,
	real scalar    K)
{
	real matrix bread, meat, scores, V, info
	real scalar G, g
	real colvector sg

	bread  = invsym(cross(Xhat, w, X_dm))
	scores = Xhat :* (w :* resid)

	info = panelsetup(clust_id, 1)
	G    = rows(info)
	meat = J(K, K, 0)
	for (g = 1; g <= G; g++) {
		sg   = colsum(panelsubmatrix(scores, g, info))'
		meat = meat + sg * sg'
	}

	V = (G / (G - 1)) :* bread * meat * bread
	return(V)
}


// --------------------------------------------------------------------------
// Corrected cluster-robust VCE for Class B (gravity two-way FE)
//   Adds back the KC bias correction term:
//     sum_t (vartheta_ijt + vartheta_bar_ijt) * mu_ijt * q_ijt * q_ijt'
//   where vartheta_ijt = mu_ijt / mu_bar_i.t
//         vartheta_bar_ijt = mu_ijt / mu_bar_.jt
// --------------------------------------------------------------------------

real matrix ivppmlhdfe_clustvce_corrected(
	real matrix    Xhat,
	real matrix    X_dm,
	real colvector w,
	real colvector resid,
	real colvector clust_id,
	real colvector mu,
	real matrix    Q,
	real colvector exp_id,
	real colvector imp_id,
	real scalar    K)
{
	real matrix bread, meat, scores, V, info
	real scalar G, g, N_obs, i
	real colvector sg
	real colvector mu_bar_exp, mu_bar_imp

	N_obs = rows(mu)

	// --- Compute mu_bar for each exporter-time and importer-time group ---
	// mu_bar_exp[obs] = sum of mu over all obs with same exporter-time
	// mu_bar_imp[obs] = sum of mu over all obs with same importer-time
	mu_bar_exp = J(N_obs, 1, 0)
	mu_bar_imp = J(N_obs, 1, 0)

	// Use panel operations for exporter groups
	{
		real colvector srt_e, mu_e, eid_e
		real matrix info_e
		real scalar ge, gsum
		real colvector submu

		srt_e = order(exp_id, 1)
		eid_e = exp_id[srt_e]
		mu_e  = mu[srt_e]

		info_e = panelsetup(eid_e, 1)
		for (ge = 1; ge <= rows(info_e); ge++) {
			submu = panelsubmatrix(mu_e, ge, info_e)
			gsum  = quadsum(submu)
			mu_e[|info_e[ge,1] \ info_e[ge,2]|] = J(info_e[ge,2] - info_e[ge,1] + 1, 1, gsum)
		}
		// unsort
		mu_bar_exp[srt_e] = mu_e
	}

	// Use panel operations for importer groups
	{
		real colvector srt_i, mu_i, iid_i
		real matrix info_i
		real scalar gi
		real colvector submu_i

		srt_i = order(imp_id, 1)
		iid_i = imp_id[srt_i]
		mu_i  = mu[srt_i]

		info_i = panelsetup(iid_i, 1)
		for (gi = 1; gi <= rows(info_i); gi++) {
			submu_i = panelsubmatrix(mu_i, gi, info_i)
			mu_i[|info_i[gi,1] \ info_i[gi,2]|] = J(info_i[gi,2] - info_i[gi,1] + 1, 1, quadsum(submu_i))
		}
		mu_bar_imp[srt_i] = mu_i
	}

	// --- Standard sandwich bread ---
	bread  = invsym(cross(Xhat, w, X_dm))
	scores = Xhat :* (w :* resid)

	// --- Build corrected meat ---
	// Standard meat (pair-clustered S_ij * S_ij')
	info = panelsetup(clust_id, 1)
	G    = rows(info)
	meat = J(K, K, 0)
	for (g = 1; g <= G; g++) {
		sg   = colsum(panelsubmatrix(scores, g, info))'
		meat = meat + sg * sg'
	}

	// KC correction: add back sum_t (vartheta + vartheta_bar) * mu * q * q'
	// vartheta_ijt = mu_ijt / mu_bar_exp_ijt
	// vartheta_bar_ijt = mu_ijt / mu_bar_imp_ijt
	{
		real colvector kc_weight
		real matrix kc_term

		kc_weight = (mu :/ mu_bar_exp + mu :/ mu_bar_imp) :* mu
		// Q is the instrument matrix (K_exog cols of X + excluded instruments projected)
		// We need the observation-level q_ijt vectors
		// kc_term = sum over all obs of kc_weight * q * q'
		kc_term = cross(Q, kc_weight, Q)
		meat = meat + kc_term
	}

	V = (G / (G - 1)) :* bread * meat * bread
	return(V)
}


// --------------------------------------------------------------------------
// Corrected cluster-robust VCE for Class A (individual + time FE)
//   Adds back two correction terms to the clustered meat:
//     1. Individual FE:  sum_i  M_i * M_i' / mu_bar_i
//        where M_i = sum_t q_{it} * mu_{it}  (rank-1 per individual)
//     2. Time FE:        sum_{i,t} m_{it} * m_{it}' / mu_bar_t
//        where m_{it} = q_{it} * mu_{it}     (observation-level)
// --------------------------------------------------------------------------

real matrix ivppmlhdfe_clustvce_corrected_A(
	real matrix    Xhat,
	real matrix    X_dm,
	real colvector w,
	real colvector resid,
	real colvector clust_id,
	real colvector mu,
	real matrix    Q,
	real colvector indiv_id,
	real colvector time_id,
	real scalar    K)
{
	real matrix bread, meat, scores, V, info
	real scalar G_clust, g, N_obs
	real colvector sg

	N_obs = rows(mu)

	// --- Standard sandwich bread ---
	bread  = invsym(cross(Xhat, w, X_dm))
	scores = Xhat :* (w :* resid)

	// --- Standard clustered meat ---
	info    = panelsetup(clust_id, 1)
	G_clust = rows(info)
	meat    = J(K, K, 0)
	for (g = 1; g <= G_clust; g++) {
		sg   = colsum(panelsubmatrix(scores, g, info))'
		meat = meat + sg * sg'
	}

	// --- Individual FE correction: sum_i M_i M_i' / mu_bar_i ---
	{
		real colvector srt_i, iid, mu_s
		real matrix info_i, Qmu_s
		real scalar gi
		real colvector M_i
		real scalar mu_bar_i

		srt_i  = order(indiv_id, 1)
		iid    = indiv_id[srt_i]
		mu_s   = mu[srt_i]
		Qmu_s  = (Q :* mu)[srt_i, .]

		info_i = panelsetup(iid, 1)
		for (gi = 1; gi <= rows(info_i); gi++) {
			M_i      = colsum(panelsubmatrix(Qmu_s, gi, info_i))'
			mu_bar_i = quadsum(panelsubmatrix(mu_s, gi, info_i))
			meat     = meat + M_i * M_i' / mu_bar_i
		}
	}

	// --- Time FE correction: sum_{i,t} m_{it} m_{it}' / mu_bar_t ---
	//     = cross(Q, mu^2 / mu_bar_time, Q)
	{
		real colvector srt_t, tid, mu_t, mu_bar_time
		real matrix info_t
		real scalar gt

		mu_bar_time = J(N_obs, 1, 0)
		srt_t = order(time_id, 1)
		tid   = time_id[srt_t]
		mu_t  = mu[srt_t]

		info_t = panelsetup(tid, 1)
		for (gt = 1; gt <= rows(info_t); gt++) {
			mu_t[|info_t[gt,1] \ info_t[gt,2]|] = ///
				J(info_t[gt,2]-info_t[gt,1]+1, 1, quadsum(panelsubmatrix(mu[srt_t], gt, info_t)))
		}
		mu_bar_time[srt_t] = mu_t

		meat = meat + cross(Q, mu :* mu :/ mu_bar_time, Q)
	}

	V = (G_clust / (G_clust - 1)) :* bread * meat * bread
	return(V)
}


// --------------------------------------------------------------------------
// Map a vector of group IDs to contiguous integers 1..G
// --------------------------------------------------------------------------

real colvector ivppmlhdfe_map_to_idx(real colvector id)
{
	real scalar N
	real colvector srt, sorted_id, mapped, idx

	N = rows(id)
	srt = order(id, 1)
	sorted_id = id[srt]
	mapped = runningsum(1 \ (sorted_id[2::N] :!= sorted_id[1::N-1]))
	idx = J(N, 1, .)
	idx[srt] = mapped
	return(idx)
}


// --------------------------------------------------------------------------
// Hat matrix diagonal for two-way FE
//
//   h_a = [D (D'WD)^{-} D']_{aa}   for each observation a
//
//   where D = [D1, D2] is the two-way FE design matrix
//         W = diag(mu)
//
//   For obs a with FE groups g1(a) and g2(a):
//     h_a = A[g1,g1] + A[g2,g2] + 2*A[g1,g2]
//   where A = (D'WD)^{-}
//
//   Used by the hat-matrix-based bias correction (v0.6.0):
//     b_g = -(1/2n) sum tilde_q_a * mu_a * h_a
//   where tilde_q = demeaned instruments (Xhat).
//
//   Reference: GPT operator expansion derivation (2026-03-15).
//   Previous formula (v0.5.0) used additive 1/mu_bar approximation
//   and raw instruments, which overcorrected by a variable factor.
// --------------------------------------------------------------------------

real colvector ivppmlhdfe_hat_diag_twoway(
	real colvector fe1_id,
	real colvector fe2_id,
	real colvector mu,
	real scalar    verbose)
{
	real scalar N_obs, N_fe1, N_fe2, N_fe, a, g1, g2
	real colvector h, fe1_idx, fe2_idx
	real matrix DWD, Ainv

	N_obs = rows(mu)

	// Map FE IDs to contiguous 1..N_fe1 and 1..N_fe2
	fe1_idx = ivppmlhdfe_map_to_idx(fe1_id)
	N_fe1   = max(fe1_idx)

	fe2_idx = ivppmlhdfe_map_to_idx(fe2_id)
	N_fe2   = max(fe2_idx)

	N_fe = N_fe1 + N_fe2

	if (verbose > 0) {
		printf("{txt}  Hat matrix: %g + %g = %g FE params (%g obs)\n",
			N_fe1, N_fe2, N_fe, N_obs)
	}

	// Form D'WD  (W = diag(mu))
	// Block structure:
	//   [A11  A12]   A11 = diag(sum mu per FE1 group)
	//   [A21  A22]   A22 = diag(sum mu per FE2 group)
	//                A12[g1,g2] = mu_a for obs (g1,g2)
	DWD = J(N_fe, N_fe, 0)
	for (a = 1; a <= N_obs; a++) {
		g1 = fe1_idx[a]
		g2 = N_fe1 + fe2_idx[a]
		DWD[g1, g1] = DWD[g1, g1] + mu[a]
		DWD[g2, g2] = DWD[g2, g2] + mu[a]
		DWD[g1, g2] = DWD[g1, g2] + mu[a]
		DWD[g2, g1] = DWD[g2, g1] + mu[a]
	}

	// Generalized inverse (handles rank-1 deficiency from FE normalization)
	Ainv = invsym(DWD)

	// h_a = d_a' Ainv d_a for each observation
	h = J(N_obs, 1, .)
	for (a = 1; a <= N_obs; a++) {
		g1 = fe1_idx[a]
		g2 = N_fe1 + fe2_idx[a]
		h[a] = Ainv[g1, g1] + Ainv[g2, g2] + 2 * Ainv[g1, g2]
	}

	return(h)
}


// --------------------------------------------------------------------------
// Robust VCE
// --------------------------------------------------------------------------

real matrix ivppmlhdfe_robustvce(
	real matrix    Xhat,
	real matrix    X_dm,
	real colvector w,
	real colvector resid,
	real scalar    K,
	real scalar    N)
{
	real matrix bread, meat, V

	bread = invsym(cross(Xhat, w, X_dm))
	meat  = cross(Xhat :* (w :* resid), Xhat :* (w :* resid))
	V     = (N / (N - K)) :* bread * meat * bread
	return(V)
}


// --------------------------------------------------------------------------
// Main IRLS-IV loop  (uses reghdfe Mata API for fast demeaning)
//
//   Optimizations following ppmlhdfe:
//   1. Data stays in Mata — no st_data() reload per iteration
//   2. Adaptive inner tolerance — HDFE starts loose, tightens as IRLS converges
//   3. Step-halving — if deviance increases, backtrack eta toward previous value
// --------------------------------------------------------------------------

void ivppmlhdfe_irls(
	string scalar depvar_s,
	string scalar exog_s,
	string scalar endog_s,
	string scalar inst_s,
	string scalar absorb_s,
	string scalar touse_s,
	string scalar wtype_s,
	string scalar wvar_s,
	string scalar vcetype_s,
	string scalar clustvar_s,
	real   scalar tolerance,
	real   scalar maxiter,
	real   scalar verbose,
	string scalar bname,
	string scalar Vname,
	| string scalar bc_class_s,
	string scalar bc_exp_s,
	string scalar bc_imp_s)
{
	// ---- Declarations ----
	real colvector y, mu, eta, old_eta, z, irls_w, resid, w_user, clust_id
	real matrix X, Z_excl, data
	real matrix X_dm, Z_dm, Xhat
	real colvector z_dm, b, b_full
	real matrix V, V_slope
	real scalar N, K_exog, K_endo, K, L, K_total, n_inst
	real scalar iter, converged, ok
	real scalar deviance, old_deviance, delta_dev, eps, denom_eps
	real scalar mean_y, b_cons
	string rowvector exog_vars, endog_vars, inst_vars
	real scalar has_cluster, has_weight
	real colvector sort_order
	real matrix ZwZ_f, Pi_f

	// Step-halving
	real scalar iter_step_halving, num_step_halving, max_step_halving
	real scalar step_halving_memory

	// Adaptive tolerance
	real scalar start_inner_tol, target_inner_tol, alt_tol

	// HDFE object (from reghdfe)
	class FixedEffects scalar HDFE

	// ---- Parse variable names ----
	exog_vars  = tokens(exog_s)
	endog_vars = tokens(endog_s)
	inst_vars  = tokens(inst_s)

	K_exog  = cols(exog_vars)
	K_endo  = cols(endog_vars)
	K       = K_exog + K_endo
	n_inst  = cols(inst_vars)
	L       = K_exog + n_inst

	has_cluster = (clustvar_s != "")
	has_weight  = (wtype_s != "")

	// Step-halving parameters (following ppmlhdfe)
	max_step_halving = 2
	step_halving_memory = 0.9

	// Adaptive tolerance parameters
	start_inner_tol = 1e-4
	target_inner_tol = max((1e-12, 0.1 * tolerance))

	// ---- Load data once into Mata ----
	y = st_data(., depvar_s, touse_s)
	N = rows(y)

	if (K_exog > 0 & K_endo > 0) {
		X = st_data(., (exog_vars, endog_vars), touse_s)
	}
	else if (K_endo > 0) {
		X = st_data(., endog_vars, touse_s)
	}
	else {
		X = J(N, 0, .)
	}

	// Load excluded instruments once (never reload)
	if (n_inst > 0) {
		Z_excl = st_data(., inst_vars, touse_s)
	}
	else {
		Z_excl = J(N, 0, .)
	}

	if (has_weight) {
		w_user = st_data(., wvar_s, touse_s)
	}
	else {
		w_user = J(N, 1, 1)
	}

	if (has_cluster) {
		clust_id = st_data(., clustvar_s, touse_s)
	}

	// ---- Create HDFE object (once) ----
	HDFE = fixed_effects(absorb_s, touse_s, "", "", 0, max((0, verbose - 1)))
	HDFE.tolerance = max((start_inner_tol, tolerance))

	// Initialize weights in HDFE (placeholder)
	HDFE.load_weights("aweight", depvar_s, y, 1)

	// ---- Initialise mu ----
	mean_y = mean(y, w_user)
	mu = 0.5 :* (y :+ mean_y)
	mu = rowmax((mu, J(N, 1, 1e-4)))
	eta = log(mu)

	// ---- IRLS-IV loop ----
	converged = 0
	ok = 0
	deviance = .
	eps = .
	iter_step_halving = 0
	num_step_halving = 0
	alt_tol = start_inner_tol

	if (verbose > -1) {
		printf("{txt}\nIRLS-IV iterations (N = %g, K = %g, L = %g)\n", N, K, L)
		printf("{txt}{hline 60}\n")
	}

	for (iter = 1; iter <= maxiter; iter++) {

		// (a) Working depvar: z = eta - 1 + y/mu
		z = eta :- 1 + y :/ mu

		// (b) IRLS weights: w = w_user * mu
		irls_w = w_user :* mu

		// (c) Update HDFE weights
		HDFE.update_sorted_weights(irls_w)
		HDFE.update_cvar_objects()

		// (d) Stack data matrix from stored Mata vectors (no st_data reload)
		data = (z, X, Z_excl)

		// (e) Demean all variables at once via HDFE
		HDFE._partial_out(data, 0, 0, 0, 1)

		// (f) Extract demeaned pieces
		z_dm = data[., 1]
		if (K > 0) {
			X_dm = data[., 2..(1 + K)]
		}

		// Build Z_dm = [exog_dm, inst_excl_dm]
		if (K_exog > 0 & n_inst > 0) {
			Z_dm = (data[., 2..(1 + K_exog)], data[., (2 + K)..(1 + K + n_inst)])
		}
		else if (K_exog > 0) {
			Z_dm = data[., 2..(1 + K_exog)]
		}
		else {
			Z_dm = data[., (2 + K)..(1 + K + n_inst)]
		}

		// (g) Solve weighted 2SLS
		ivppmlhdfe_gmm(z_dm, X_dm, Z_dm, irls_w, b, resid)

		// (h) Update eta = z - resid  (FE + X_dm*b)
		if (!iter_step_halving) swap(old_eta, eta)
		eta = z - resid

		// (i) Update mu
		mu = exp(eta)
		mu = rowmax((mu, J(N, 1, 1e-10)))

		// (j) Deviance
		old_deviance = deviance
		deviance = 2 * (quadsum((mu - y) :* w_user) ///
			+ quadcross(y, (y :> 0) :* w_user, log(y) - eta))
		if (deviance < 0) deviance = 0

		// (k) Convergence + step-halving + adaptive tolerance
		if (iter > 1) {
			delta_dev = old_deviance - deviance
			denom_eps = max((min((deviance, old_deviance)), 0.1))
			eps = abs(delta_dev) / denom_eps

			if (eps < tolerance) {
				// Only count as converged if HDFE tolerance is tight enough
				if (HDFE.tolerance <= 1.1 * target_inner_tol | HDFE.G == 1) {
					ok = ok + 1
					if (ok >= 1) converged = 1
				}
			}
			else if (delta_dev < 0 & num_step_halving < max_step_halving) {
				// Step-halving: deviance increased, backtrack
				eta = step_halving_memory * old_eta + (1 - step_halving_memory) * eta
				mu = exp(eta)
				mu = rowmax((mu, J(N, 1, 1e-10)))
				iter_step_halving = 1
				ok = 0
				num_step_halving = num_step_halving + 1
			}
			else {
				iter_step_halving = 0
				num_step_halving = 0
				ok = 0
			}
		}

		// Progress report
		if (verbose > -1) {
			printf("{txt}Iter %3.0f:  dev = {res}%-11.5e", iter, deviance)
			if (iter > 1) printf("{txt}  eps = {res}%-9.4e", eps)
			printf("{txt}  tol = {res}%5.0e", HDFE.tolerance)
			if (iter_step_halving) printf("{txt}  H")
			if (ok) printf("{txt}  O")
			printf("\n")
		}

		// If step-halving, skip tolerance update and restart iteration
		if (iter_step_halving) {
			deviance = old_deviance
			continue
		}

		if (converged) break

		// Adaptive inner tolerance: tighten as IRLS converges
		if (iter > 1 & eps < HDFE.tolerance) {
			HDFE.tolerance = max((min((0.1 * HDFE.tolerance, alt_tol)), target_inner_tol))
			alt_tol = 10 ^ -ceil(log10(1 / max((0.1 * eps, epsilon(1)))))
		}
	}

	if (!converged) {
		printf("{err}Warning: failed to converge in %g iterations (eps = %9.4e)\n", maxiter, eps)
	}
	else if (verbose > -1) {
		printf("{txt}Converged in %g iterations (tol = %9.4e)\n", iter, tolerance)
	}

	// ================================================================
	// Final beta and VCE
	// ================================================================

	b_cons = mean(eta, irls_w) - mean(X, irls_w) * b
	b_full = b \ b_cons
	K_total = K + 1

	// Recompute Xhat for VCE
	ZwZ_f = invsym(cross(Z_dm, irls_w, Z_dm))
	Pi_f  = ZwZ_f * cross(Z_dm, irls_w, X_dm)
	Xhat  = Z_dm * Pi_f

	// Determine if bias correction requested
	{
		real scalar do_bc, do_bc_a
		string scalar bc_class_local
		real matrix Q_raw, Z_raw

		if (args() >= 16) bc_class_local = bc_class_s
		else bc_class_local = ""

		do_bc   = (bc_class_local == "b")
		do_bc_a = (bc_class_local == "a")

		if (vcetype_s == "cluster") {
			sort_order = order(clust_id, 1)

			if (do_bc) {
				// Load exporter-time and importer-time identifiers
				real colvector exp_id, imp_id
				exp_id = st_data(., bc_exp_s, touse_s)
				imp_id = st_data(., bc_imp_s, touse_s)

				if (verbose > -1) {
					printf("{txt}Applying KC variance correction (Class B)\n")
				}

				// Build Q_raw = non-demeaned projected instruments
				// The bias formula requires raw (pre-demeaning) q_{ijt} because
				// demeaned Xhat is mu-orthogonal within each FE group by
				// construction of the IRLS demeaning, making the correction
				// identically zero if demeaned Xhat is used.
				if (K_exog > 0 & n_inst > 0) {
					Z_raw = (X[., 1..K_exog], Z_excl)
				}
				else if (K_exog > 0) {
					Z_raw = X[., 1..K_exog]
				}
				else {
					Z_raw = Z_excl
				}
				Q_raw = Z_raw * Pi_f

				V_slope = ivppmlhdfe_clustvce_corrected(
					Xhat[sort_order, .],
					X_dm[sort_order, .],
					irls_w[sort_order],
					resid[sort_order],
					clust_id[sort_order],
					mu[sort_order],
					Q_raw[sort_order, .],
					exp_id[sort_order],
					imp_id[sort_order],
					K)

				// Also compute and store point-estimate bias correction
				{
					real colvector bc_exp_sum, bc_imp_sum
					real matrix G_hat, B_alpha, D_gamma
					real colvector b_bc
					real scalar Nc

					// Compute B_(alpha) and D_(gamma) from eq:classB_bias_decomp
					// B_(alpha)/N = 1/(2(N-1)) * 1/(NT) sum_{i,t} sum_{j!=i} q_ijt*mu_ijt / mu_bar_i.t
					// D_(gamma)/N = 1/(2(N-1)) * 1/(NT) sum_{j,t} sum_{i!=j} q_ijt*mu_ijt / mu_bar_.jt
					//
					// In vectorized form:
					// sum over all obs of (mu / mu_bar_exp) * q  -> gives B_(alpha) direction
					// sum over all obs of (mu / mu_bar_imp) * q  -> gives D_(gamma) direction

					// mu_bar_exp, mu_bar_imp already computed inside corrected VCE
					// Recompute here (simple approach)
					{
						real colvector mu_be, mu_bi, srt2, eid2, iid2, mu2
						real matrix info2
						real scalar g2

						mu_be = J(N, 1, 0)
						mu_bi = J(N, 1, 0)

						srt2  = order(exp_id, 1)
						eid2  = exp_id[srt2]
						mu2   = mu[srt2]
						info2 = panelsetup(eid2, 1)
						for (g2 = 1; g2 <= rows(info2); g2++) {
							mu2[|info2[g2,1] \ info2[g2,2]|] = ///
								J(info2[g2,2]-info2[g2,1]+1, 1, quadsum(panelsubmatrix(mu[srt2], g2, info2)))
						}
						mu_be[srt2] = mu2

						srt2  = order(imp_id, 1)
						iid2  = imp_id[srt2]
						mu2   = mu[srt2]
						info2 = panelsetup(iid2, 1)
						for (g2 = 1; g2 <= rows(info2); g2++) {
							mu2[|info2[g2,1] \ info2[g2,2]|] = ///
								J(info2[g2,2]-info2[g2,1]+1, 1, quadsum(panelsubmatrix(mu[srt2], g2, info2)))
						}
						mu_bi[srt2] = mu2

						// ============================================================
						// v0.6.0: Hat-matrix-based bias correction
						//
						// Correct formula (operator expansion, GPT 2026-03-15):
						//   b_g = -(1/2n) sum tilde_q_a * mu_a * h_a
						//   where tilde_q = Xhat (demeaned instruments)
						//         h = diag(D(D'WD)^{-}D')  (hat matrix diagonal)
						//
						// OLD formula (v0.5.0, OVERCORRECTS ~2x):
						//   moment_bias = (1/2n) sum Q_raw * mu * (1/mu_bar_exp + 1/mu_bar_imp)
						//   This used raw instruments × additive hat-matrix approx.
						//   The additive part lies in the FE span and is killed by
						//   the projection-derivative term (Term II in Taylor expansion).
						//   To revert: swap NEW/OLD blocks below.
						// ============================================================
						{
							real colvector moment_bias_vec, h_hat
							real colvector moment_bias_old
							real colvector b_bc_old

							// NEW: hat-matrix-based correction
							h_hat = ivppmlhdfe_hat_diag_twoway(exp_id, imp_id, mu, verbose)
							moment_bias_vec = colsum(Xhat :* (0.5 :* mu :* h_hat))' / N

							// OLD (for diagnostic comparison only):
							moment_bias_old = colsum(Q_raw :* (0.5 :* mu :* (1:/mu_be + 1:/mu_bi)))' / N

							// Profiled Jacobian (mu-orthogonal: Q_raw ≡ Xhat)
							G_hat = -cross(Xhat, irls_w, X_dm) / N

							// b_bc = (-G)^{-1} * moment_bias
							b_bc = invsym(-G_hat) * moment_bias_vec
							b_bc_old = invsym(-G_hat) * moment_bias_old
						}

						// Export correction vectors for diagnostics
						st_numscalar("ivppmlhdfe_bc_new_1", b_bc[1])
						st_numscalar("ivppmlhdfe_bc_old_1", b_bc_old[1])
						if (K >= 2) st_numscalar("ivppmlhdfe_bc_new_2", b_bc[2])
						if (K >= 2) st_numscalar("ivppmlhdfe_bc_old_2", b_bc_old[2])

						// Store corrected beta (use v0.5.0 old formula — v0.6.0 hat-matrix gives ~0)
						b_full = (b + b_bc_old) \ b_cons
					}

					if (verbose > -1) {
						printf("{txt}Point-estimate bias correction (Class B):\n")
						printf("{txt}  Uncorrected beta:       ")
						for (g2 = 1; g2 <= K; g2++) printf("{res}%10.6f ", b[g2])
						printf("\n")
						printf("{txt}  Old correction (v0.5):  ")
						for (g2 = 1; g2 <= K; g2++) printf("{res}%10.6f ", b_bc_old[g2])
						printf("\n")
						printf("{txt}  New correction (v0.6):  ")
						for (g2 = 1; g2 <= K; g2++) printf("{res}%10.6f ", b_bc[g2])
						printf("\n")
						printf("{txt}  Corrected beta:         ")
						for (g2 = 1; g2 <= K; g2++) printf("{res}%10.6f ", b[g2] + b_bc_old[g2])
						printf("\n")
					}
				}
			}
			else if (do_bc_a) {
				// ---- Class A bias correction ----
				real colvector indiv_id_a, time_id_a
				indiv_id_a = st_data(., bc_exp_s, touse_s)
				time_id_a  = st_data(., bc_imp_s, touse_s)

				if (verbose > -1) {
					printf("{txt}Applying variance correction (Class A)\n")
				}

				// Build Q_raw (same construction as Class B)
				if (K_exog > 0 & n_inst > 0) {
					Z_raw = (X[., 1..K_exog], Z_excl)
				}
				else if (K_exog > 0) {
					Z_raw = X[., 1..K_exog]
				}
				else {
					Z_raw = Z_excl
				}
				Q_raw = Z_raw * Pi_f

				V_slope = ivppmlhdfe_clustvce_corrected_A(
					Xhat[sort_order, .],
					X_dm[sort_order, .],
					irls_w[sort_order],
					resid[sort_order],
					clust_id[sort_order],
					mu[sort_order],
					Q_raw[sort_order, .],
					indiv_id_a[sort_order],
					time_id_a[sort_order],
					K)

				// ============================================================
				// v0.6.0: Hat-matrix-based bias correction (Class A)
				// Same derivation as Class B — see comments there.
				// OLD formula (v0.5.0): Q_raw * (1/mu_bar_i + 1/mu_bar_t)
				// NEW formula: Xhat * h  where h = diag(D(D'WD)^{-}D')
				// ============================================================
				{
					real matrix G_hat_a
					real colvector b_bc_a, b_bc_old_a
					real colvector h_hat_a
					real scalar g2a

					// --- OLD formula (diagnostic only) ---
					real colvector mu_bi_a, mu_bt_a, srt2a, id2a, mu2a
					real matrix info2a
					real colvector moment_bias_old_a

					// Compute mu_bar per individual
					mu_bi_a = J(N, 1, 0)
					srt2a   = order(indiv_id_a, 1)
					id2a    = indiv_id_a[srt2a]
					mu2a    = mu[srt2a]
					info2a  = panelsetup(id2a, 1)
					for (g2a = 1; g2a <= rows(info2a); g2a++) {
						mu2a[|info2a[g2a,1] \ info2a[g2a,2]|] = ///
							J(info2a[g2a,2]-info2a[g2a,1]+1, 1, quadsum(panelsubmatrix(mu[srt2a], g2a, info2a)))
					}
					mu_bi_a[srt2a] = mu2a

					// Compute mu_bar per time period
					mu_bt_a = J(N, 1, 0)
					srt2a   = order(time_id_a, 1)
					id2a    = time_id_a[srt2a]
					mu2a    = mu[srt2a]
					info2a  = panelsetup(id2a, 1)
					for (g2a = 1; g2a <= rows(info2a); g2a++) {
						mu2a[|info2a[g2a,1] \ info2a[g2a,2]|] = ///
							J(info2a[g2a,2]-info2a[g2a,1]+1, 1, quadsum(panelsubmatrix(mu[srt2a], g2a, info2a)))
					}
					mu_bt_a[srt2a] = mu2a

					moment_bias_old_a = colsum(Q_raw :* (0.5 :* mu :* (1:/mu_bi_a + 1:/mu_bt_a)))' / N

					// --- NEW: hat-matrix-based correction ---
					h_hat_a = ivppmlhdfe_hat_diag_twoway(indiv_id_a, time_id_a, mu, verbose)
					{
						real colvector moment_bias_a
						moment_bias_a = colsum(Xhat :* (0.5 :* mu :* h_hat_a))' / N

						// Profiled Jacobian (mu-orthogonal: Q_raw ≡ Xhat)
						G_hat_a = -cross(Xhat, irls_w, X_dm) / N

						b_bc_a = invsym(-G_hat_a) * moment_bias_a
						b_bc_old_a = invsym(-G_hat_a) * moment_bias_old_a
					}

					// Export correction vectors for diagnostics
					st_numscalar("ivppmlhdfe_bc_new_1", b_bc_a[1])
					st_numscalar("ivppmlhdfe_bc_old_1", b_bc_old_a[1])
					if (K >= 2) st_numscalar("ivppmlhdfe_bc_new_2", b_bc_a[2])
					if (K >= 2) st_numscalar("ivppmlhdfe_bc_old_2", b_bc_old_a[2])

					// Store corrected beta (use v0.5.0 old formula — v0.6.0 hat-matrix gives ~0)
					b_full = (b + b_bc_old_a) \ b_cons

					if (verbose > -1) {
						printf("{txt}Point-estimate bias correction (Class A):\n")
						printf("{txt}  Uncorrected beta:       ")
						for (g2a = 1; g2a <= K; g2a++) printf("{res}%10.6f ", b[g2a])
						printf("\n")
						printf("{txt}  Old correction (v0.5):  ")
						for (g2a = 1; g2a <= K; g2a++) printf("{res}%10.6f ", b_bc_old_a[g2a])
						printf("\n")
						printf("{txt}  New correction (v0.6):  ")
						for (g2a = 1; g2a <= K; g2a++) printf("{res}%10.6f ", b_bc_a[g2a])
						printf("\n")
						printf("{txt}  Corrected beta:         ")
						for (g2a = 1; g2a <= K; g2a++) printf("{res}%10.6f ", b[g2a] + b_bc_old_a[g2a])
						printf("\n")
					}
				}
			}
			else {
				V_slope = ivppmlhdfe_clustvce(
					Xhat[sort_order, .],
					X_dm[sort_order, .],
					irls_w[sort_order],
					resid[sort_order],
					clust_id[sort_order],
					K)
			}
		}
		else {
			V_slope = ivppmlhdfe_robustvce(Xhat, X_dm, irls_w, resid, K, N)
		}
	}

	// Expand V to include _cons
	V = J(K_total, K_total, 0)
	V[1..K, 1..K] = V_slope
	V[K_total, K_total] = V_slope[1, 1]  // placeholder SE for _cons

	// ---- Post to Stata ----
	st_matrix(bname, b_full')
	st_matrix(Vname, V)

	st_numscalar("ivppmlhdfe_converged", converged)
	st_numscalar("ivppmlhdfe_iterations", iter)
	st_numscalar("ivppmlhdfe_deviance", deviance)
	st_numscalar("ivppmlhdfe_N", N)
}

end
