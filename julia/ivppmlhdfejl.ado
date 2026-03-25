*! ivppmlhdfejl 0.8.0  24mar2026
*! IV-PPML with High-Dimensional Fixed Effects — Julia backend
*! Architecture: Path 3 — single Julia call with solver reuse (update_weights!)
*! Mirrors reghdfejl / GLFixedEffectModels.jl pattern
*!
*! Algorithm: IRLS-IV  (Mullahy 1997 GMM via iteratively reweighted 2SLS)
*! FE absorption via FixedEffects.jl (iterative demeaning)
*! Authors: Ohyun Kwon
*!
*! Syntax:
*!   ivppmlhdfejl depvar [exogvars] (endogvars = excluded_instruments) [if] [in] [pw], ///
*!       absorb(absvars) [vce(cluster c1 [c2]) tolerance(#) maxiterations(#) verbose(#) ///
*!       exposure(varname) offset(varname) d(name) separation(string)]

program ivppmlhdfejl, eclass
	version 15.0

	if replay() {
		if "`e(cmd)'" != "ivppmlhdfejl" error 301
		_coef_table, `0'
		exit
	}

	* Initialize Julia environment
	qui jl GetEnv
	local env `r(env)'
	ivppmlhdfejl_load

	cap noi _ivppmlhdfejl `0'
	local rc = _rc

	* Restore previous Julia environment
	qui jl SetEnv `env'

	exit `rc'
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
		di as err "Syntax: ivppmlhdfejl depvar [exogvars] (endogvars = instruments) [if] [in], absorb() [options]"
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

program _ivppmlhdfejl, eclass
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
		 EXPosure(varname) OFFset(varname) ///
		 D(name) SEParation(string) ///
		 gpu THReads(integer 0) KEEPSINgletons compact]

	if "`irr'" != "" local eform "eform"

	// exposure and offset
	if "`exposure'" != "" & "`offset'" != "" {
		di as err "cannot specify both exposure() and offset()"
		exit 198
	}
	local offset_display ""
	if "`exposure'" != "" {
		tempvar offset_var
		quietly gen double `offset_var' = ln(`exposure')
		local offset `offset_var'
		local offset_display "ln(`exposure')"
	}
	else if "`offset'" != "" {
		local offset_display "`offset'"
	}

	// d() option
	if "`d'" != "" {
		confirm new variable `d'
	}

	// separation
	if inlist("`separation'", "", "def", "default", "on", "auto") {
		local separation "fe"
	}
	if inlist("`separation'", "no", "off", "none") {
		local separation ""
	}
	local do_separation = ("`separation'" != "")

	// cluster / vce parsing
	if "`cluster'" != "" {
		if "`vce'" != "" {
			di as err "cannot specify both cluster() and vce()"
			exit 198
		}
		local vce "cluster `cluster'"
	}

	local vcetype "robust"
	local clustvars ""
	local n_clustvars = 0
	if "`vce'" != "" {
		gettoken vtype vrest : vce
		if "`vtype'" == "cluster" {
			local vcetype "cluster"
			local clustvars = strtrim(`"`vrest'"')
			local n_clustvars : word count `clustvars'
			if `n_clustvars' == 0 {
				di as err "vce(cluster) requires at least one variable"
				exit 198
			}
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
	if "`clustvars'" != "" {
		local clustvars_unab ""
		foreach cv of local clustvars {
			unab cv : `cv'
			local clustvars_unab `clustvars_unab' `cv'
		}
		local clustvars `clustvars_unab'
	}

	local regressors   `exog' `endog'
	local n_exog     : word count `exog'
	local n_endog    : word count `endog'
	local n_inst     : word count `instruments'

	if `n_inst' < `n_endog' {
		di as err "equation not identified: need at least `n_endog' instruments, have `n_inst'"
		exit 481
	}

	// mark sample
	tempvar touse
	mark `touse' `if' `in'
	markout `touse' `depvar' `regressors' `instruments' `clustvars'
	if "`weight'" != "" markout `touse' `exp'
	if "`offset'" != "" markout `touse' `offset'

	quietly count if `depvar' < 0 & `touse'
	if r(N) > 0 {
		di as err "`depvar' must be nonnegative"
		exit 459
	}

	quietly count if `touse'
	local N_full = r(N)
	if `N_full' == 0 error 2000

	// ---------- 4.  Separation detection (Stata-side, FE method) ----------
	local n_separated = 0
	if `do_separation' {
		local sep_changed = 1
		while `sep_changed' {
			local sep_changed = 0
			local absorb_clean : subinstr local absorb "#" " ", all
			foreach tok of local absorb_clean {
				local v = regexr("`tok'", "^[ic]*\.", "")
				capture confirm numeric variable `v'
				if _rc continue
				tempvar grpsum
				quietly egen double `grpsum' = total(`depvar') if `touse', by(`v')
				quietly count if `grpsum' == 0 & `touse'
				if r(N) > 0 {
					quietly replace `touse' = 0 if `grpsum' == 0 & `touse'
					local n_separated = `n_separated' + r(N)
					local sep_changed = 1
				}
				drop `grpsum'
			}
		}
	}

	quietly count if `touse'
	local N_after_sep = r(N)
	if `N_after_sep' == 0 {
		di as err "no observations remaining after separation detection"
		exit 2000
	}

	// ---------- 5.  d() variable ----------
	if "`d'" != "" {
		quietly gen double `d' = .
	}

	// ---------- 6.  Parse absorb ----------
	local absorb_raw `"`absorb'"'

	local fe_vars ""
	local fe_syms ""
	local n_fe 0
	foreach avar of local absorb {
		local avar_clean = regexr("`avar'", "^i\.", "")
		local ++n_fe
		tempvar fe_`n_fe'
		qui egen long `fe_`n_fe'' = group(`avar_clean') if `touse'
		local fe_vars `fe_vars' `fe_`n_fe''
		local fe_syms `fe_syms' `fe_`n_fe''
	}

	// ---------- 7.  User weights ----------
	local has_weight = ("`weight'" != "")
	if `has_weight' {
		cap confirm var `exp'
		if _rc {
			tempvar wtvar
			qui gen double `wtvar' = `exp' if `touse'
		}
		else local wtvar `exp'
	}

	// ---------- 8.  Cluster variables ----------
	local has_cluster = ("`clustvars'" != "")
	if `has_cluster' {
		local clustvar_nums ""
		local j = 0
		foreach cv of local clustvars {
			local ++j
			cap confirm numeric var `cv'
			if _rc {
				tempvar clustvar_num_`j'
				qui egen long `clustvar_num_`j'' = group(`cv') if `touse'
				local clustvar_nums `clustvar_nums' `clustvar_num_`j''
			}
			else {
				local clustvar_nums `clustvar_nums' `cv'
			}
		}
	}

	// ---------- 9.  Display pre-estimation ----------
	if `verbose' > 0 {
		di as txt _n "ivppmlhdfejl: IV-PPML with HDFE (Julia)"
		di as txt "  Dep var:      `depvar'"
		di as txt "  Exogenous:    `exog'"
		di as txt "  Endogenous:   `endog'"
		di as txt "  Instruments:  `instruments'"
		di as txt "  Absorb:       `absorb'"
		di as txt "  VCE:          `vcetype' `clustvars'"
		if "`offset_display'" != "" di as txt "  Offset:       `offset_display'"
		if `n_separated' > 0 di as txt "  Separated:    `n_separated' obs dropped"
	}

	// ---------- 10.  Transfer data to Julia DataFrame ----------
	local putvars `depvar' `exog' `endog' `instruments' `fe_vars'
	if `has_weight' local putvars `putvars' `wtvar'
	if `has_cluster' local putvars `putvars' `clustvar_nums'
	if "`offset'" != "" local putvars `putvars' `offset'
	local putvars: list uniq putvars

	jl PutVarsToDF `putvars' if `touse', nomissing doubleonly nolabel

	// ---------- 11.  Single Julia call (Path 3 architecture) ----------

	// Build Julia Symbol vectors using Symbol("name") syntax
	local jl_exog ""
	if `n_exog' > 0 {
		foreach v of local exog {
			if "`jl_exog'" == "" local jl_exog Symbol("`v'")
			else                local jl_exog `jl_exog', Symbol("`v'")
		}
	}
	local jl_endog ""
	foreach v of local endog {
		if "`jl_endog'" == "" local jl_endog Symbol("`v'")
		else                 local jl_endog `jl_endog', Symbol("`v'")
	}
	local jl_inst ""
	foreach v of local instruments {
		if "`jl_inst'" == "" local jl_inst Symbol("`v'")
		else                local jl_inst `jl_inst', Symbol("`v'")
	}
	local jl_fes ""
	foreach v of local fe_syms {
		if "`jl_fes'" == "" local jl_fes Symbol("`v'")
		else               local jl_fes `jl_fes', Symbol("`v'")
	}

	// Cluster keyword argument (supports multi-way)
	local jl_clustopt ""
	if `has_cluster' {
		if `n_clustvars' == 1 {
			local cv1 : word 1 of `clustvar_nums'
			local jl_clustopt , cluster = Symbol("`cv1'")
		}
		else {
			local jl_clustsyms ""
			foreach cv of local clustvar_nums {
				if "`jl_clustsyms'" == "" local jl_clustsyms Symbol("`cv'")
				else                     local jl_clustsyms `jl_clustsyms', Symbol("`cv'")
			}
			local jl_clustopt , cluster = Symbol[`jl_clustsyms']
		}
	}

	// Weight keyword argument
	local jl_wtopt ""
	if `has_weight' {
		local jl_wtopt , weights = Symbol("`wtvar'")
	}

	// Offset keyword argument
	local jl_offopt ""
	if "`offset'" != "" {
		local jl_offopt , offset = Symbol("`offset'")
	}

	// Estimate!
	_jl: m = IVPPMLFixedEffectModels.ivppml_reg(df, ///
		Symbol("`depvar'"), ///
		Symbol[`jl_exog'], ///
		Symbol[`jl_endog'], ///
		Symbol[`jl_inst'], ///
		Symbol[`jl_fes'], ///
		tol = `tolerance', maxiter = `maxiterations', verbose = `verbose' ///
		`jl_clustopt' `jl_wtopt' `jl_offopt');

	// Free DataFrame
	_jl: df = nothing;

	// ---------- 12.  Retrieve results ----------
	tempname b_out V_out t

	_jl: st_numscalar("`t'", m.N);
	local N_final = `t'

	_jl: st_numscalar("`t'", m.converged ? 1.0 : 0.0);
	local converged_final = `t'

	_jl: st_numscalar("`t'", Float64(m.iterations));
	local iter_final = `t'

	_jl: st_numscalar("`t'", m.deviance);
	local dev_final = `t'

	_jl: st_numscalar("`t'", m.ll);
	local ll_final = `t'

	_jl: st_numscalar("`t'", m.ll_0);
	local ll0_final = `t'

	// Get coefficient vector
	_jl: _ivppml_b = collect(IVPPMLFixedEffectModels.coef(m)');
	jl GetMatFromMat `b_out', source(_ivppml_b)

	// Get VCE matrix
	_jl: _ivppml_V = IVPPMLFixedEffectModels.vcov(m);
	jl GetMatFromMat `V_out', source(_ivppml_V)

	// Get cluster counts
	if `has_cluster' {
		forvalues j = 1/`n_clustvars' {
			_jl: st_numscalar("`t'", Float64(m.N_clust[`j']));
			local G`j' = `t'
		}
	}

	// d() variable: retrieve FE sum from Julia
	if "`d'" != "" {
		_jl: _ivppml_d = reshape(m.d_values, 1, length(m.d_values));
		tempname d_mat
		jl GetMatFromMat `d_mat', source(_ivppml_d)
		mata: st_store(., st_varindex("`d'"), "`touse'", st_matrix("`d_mat'")')
		_jl: _ivppml_d = nothing;
	}

	// Clean up Julia memory
	_jl: _ivppml_b = _ivppml_V = nothing;
	_jl: m = nothing;

	// ---------- 13.  Post results ----------
	local colnames `regressors' _cons
	matrix colnames `b_out' = `colnames'
	matrix colnames `V_out' = `colnames'
	matrix rownames `V_out' = `colnames'

	ereturn post `b_out' `V_out', esample(`touse') depname(`depvar')

	// Scalars
	ereturn scalar N           = `N_final'
	ereturn scalar converged   = `converged_final'
	ereturn scalar iterations  = `iter_final'
	ereturn scalar deviance    = `dev_final'
	ereturn scalar ll          = `ll_final'
	ereturn scalar ll_0        = `ll0_final'
	ereturn scalar r2_p        = 1 - `ll_final' / `ll0_final'
	ereturn scalar N_full      = `N_full'
	ereturn scalar num_separated = `n_separated'

	local K = `n_exog' + `n_endog'
	ereturn scalar df_m        = `K'
	ereturn scalar rank        = `K' + 1

	// Wald chi2
	if `K' > 0 {
		tempname b_slope V_slope_test chi2_mat
		matrix `b_slope' = e(b)
		matrix `b_slope' = `b_slope'[1, 1..`K']
		matrix `V_slope_test' = e(V)
		matrix `V_slope_test' = `V_slope_test'[1..`K', 1..`K']
		capture matrix `V_slope_test' = syminv(`V_slope_test')
		if !_rc {
			matrix `chi2_mat' = `b_slope' * `V_slope_test' * `b_slope''
			ereturn scalar chi2 = `chi2_mat'[1,1]
		}
		else {
			ereturn scalar chi2 = .
		}
	}
	else {
		ereturn scalar chi2 = 0
	}

	// Cluster counts
	if "`vcetype'" == "cluster" {
		ereturn scalar N_clustervars = `n_clustvars'
		forvalues j = 1/`n_clustvars' {
			local cvar : word `j' of `clustvars'
			ereturn local clustvar`j' "`cvar'"
			ereturn scalar N_clust`j' = `G`j''
		}
		// e(df) = min clusters - 1
		ereturn scalar df = `G1' - 1
		if `n_clustvars' > 1 {
			local min_G = `G1'
			forvalues j = 2/`n_clustvars' {
				if `G`j'' < `min_G' local min_G = `G`j''
			}
			ereturn scalar df = `min_G' - 1
		}
		ereturn scalar N_clust = e(df) + 1
	}
	else {
		ereturn scalar df = `N_final' - `K' - 1
	}

	// Locals
	ereturn local cmd          "ivppmlhdfejl"
	ereturn local cmdline      "ivppmlhdfejl `0'"
	ereturn local depvar       "`depvar'"
	ereturn local exogvars     "`exog'"
	ereturn local endogvars    "`endog'"
	ereturn local instruments  "`instruments'"
	ereturn local absorb       "`absorb_raw'"
	ereturn local absvars      "`absorb_raw'"
	ereturn local vcetype      "`vcetype'"
	ereturn local clustvar     "`clustvars'"
	ereturn local title        "IV-PPML with HDFE (Julia)"
	ereturn local chi2type     "Wald"
	ereturn local predict      "ivppmlhdfe_p"
	ereturn local marginsok    "default"
	ereturn local properties   "b V"
	ereturn local separation   "`separation'"
	if "`offset_display'" != "" ereturn local offset "`offset_display'"
	if "`d'" != "" ereturn local d "`d'"

	// ---------- 14.  Display ----------
	di as txt _n "{hline 78}"
	di as txt "IV-PPML with High-Dimensional Fixed Effects  (Julia backend v0.8)"
	di as txt "{hline 78}"

	* Left-right header (following ppmlhdfe_header / ivppmlhdfe)
	di as txt "Dependent variable: " _col(24) as res "`depvar'" ///
		_col(51) as txt "No. of obs" _col(68) "=" _col(70) as res %10.0gc e(N)
	di as txt "Endogenous:         " _col(24) as res "`endog'" ///
		_col(51) as txt "Residual df" _col(68) "=" _col(70) as res %10.0gc e(df)
	di as txt "Instruments:        " _col(24) as res "`instruments'" ///
		_col(51) as txt "Wald chi2(`K')" _col(68) "=" _col(70) as res %10.2f e(chi2)
	di as txt "Absorbed FE:        " _col(24) as res "`absorb_raw'" ///
		_col(51) as txt "Prob > chi2" _col(68) "=" _col(70) as res %10.4f chi2tail(`K', e(chi2))

	if "`offset_display'" != "" {
		di as txt "Offset:             " _col(24) as res "`offset_display'" ///
			_col(51) as txt "Pseudo R2" _col(68) "=" _col(70) as res %10.4f e(r2_p)
	}
	else {
		di as txt _col(51) "Pseudo R2" _col(68) "=" _col(70) as res %10.4f e(r2_p)
	}

	di as txt "Log pseudolikelihood = " as res %12.0g e(ll) ///
		_col(51) as txt "Deviance" _col(68) "=" _col(70) as res %10.4g e(deviance)

	if "`vcetype'" == "cluster" {
		forvalues j = 1/`n_clustvars' {
			local cvar : word `j' of `clustvars'
			di _col(51) as txt "No. of clusters" _col(68) "=" _col(70) as res %10.0fc e(N_clust`j') ///
				as txt "  (`cvar')"
		}
	}

	if e(converged) == 1 {
		di as txt "Converged:          " as res "yes" ///
			as txt "  (iterations = " as res e(iterations) as txt ")"
	}
	else {
		di as txt "Converged:          " as err "no" ///
			as txt "  (iterations = " as res e(iterations) as txt ")"
	}

	if e(num_separated) > 0 {
		di as txt "(`=e(num_separated)' separated observations dropped)"
	}

	di as txt "{hline 78}"
	_coef_table, `eform'
	di as txt "{hline 78}"
	di as txt "Endogenous: " as res "`endog'"
	di as txt "Instruments:" as res " `instruments'"
end
