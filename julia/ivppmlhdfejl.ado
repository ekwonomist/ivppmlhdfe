*! ivppmlhdfejl 0.2.0  14mar2026
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
*!       absorb(absvars) [vce(cluster clustvar) tolerance(#) maxiterations(#) verbose(#)]

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
		 gpu THReads(integer 0) KEEPSINgletons compact]

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
	markout `touse' `depvar' `regressors' `instruments' `clustvar'
	if "`weight'" != "" markout `touse' `exp'

	quietly count if `depvar' < 0 & `touse'
	if r(N) > 0 {
		di as err "`depvar' must be nonnegative"
		exit 459
	}

	quietly count if `touse'
	if r(N) == 0 error 2000

	// ---------- 4.  Parse absorb ----------
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

	// ---------- 5.  User weights ----------
	local has_weight = ("`weight'" != "")
	if `has_weight' {
		cap confirm var `exp'
		if _rc {
			tempvar wtvar
			qui gen double `wtvar' = `exp' if `touse'
		}
		else local wtvar `exp'
	}

	// ---------- 6.  Cluster variable ----------
	local has_cluster = ("`clustvar'" != "")
	if `has_cluster' {
		cap confirm numeric var `clustvar'
		if _rc {
			tempvar clustvar_num
			qui egen long `clustvar_num' = group(`clustvar') if `touse'
		}
		else local clustvar_num `clustvar'
	}

	// ---------- 7.  Transfer data to Julia DataFrame ----------
	local putvars `depvar' `exog' `endog' `instruments' `fe_vars'
	if `has_weight' local putvars `putvars' `wtvar'
	if `has_cluster' local putvars `putvars' `clustvar_num'
	local putvars: list uniq putvars

	jl PutVarsToDF `putvars' if `touse', nomissing doubleonly nolabel

	// ---------- 8.  Single Julia call (Path 3 architecture) ----------

	// Build Julia Symbol vectors using Symbol("name") syntax
	// (avoid :name — Stata's _jl: parser interprets colon as Stata syntax)
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

	// Cluster keyword argument
	local jl_clustopt ""
	if `has_cluster' {
		local jl_clustopt , cluster = Symbol("`clustvar_num'")
	}

	// Weight keyword argument
	local jl_wtopt ""
	if `has_weight' {
		local jl_wtopt , weights = Symbol("`wtvar'")
	}

	// Estimate!
	_jl: m = IVPPMLFixedEffectModels.ivppml_reg(df, ///
		Symbol("`depvar'"), ///
		Symbol[`jl_exog'], ///
		Symbol[`jl_endog'], ///
		Symbol[`jl_inst'], ///
		Symbol[`jl_fes'], ///
		tol = `tolerance', maxiter = `maxiterations', verbose = `verbose' ///
		`jl_clustopt' `jl_wtopt');

	// Free DataFrame
	_jl: df = nothing;

	// ---------- 9.  Retrieve results ----------
	tempname b_out V_out t

	_jl: st_numscalar("`t'", m.N);
	local N_final = `t'

	_jl: st_numscalar("`t'", m.converged ? 1.0 : 0.0);
	local converged_final = `t'

	_jl: st_numscalar("`t'", Float64(m.iterations));
	local iter_final = `t'

	_jl: st_numscalar("`t'", m.deviance);
	local dev_final = `t'

	// Get coefficient vector
	_jl: _ivppml_b = collect(IVPPMLFixedEffectModels.coef(m)');
	jl GetMatFromMat `b_out', source(_ivppml_b)

	// Get VCE matrix
	_jl: _ivppml_V = IVPPMLFixedEffectModels.vcov(m);
	jl GetMatFromMat `V_out', source(_ivppml_V)

	// Clean up Julia memory
	_jl: _ivppml_b = _ivppml_V = nothing;
	_jl: m = nothing;

	// ---------- 10.  Post results ----------
	local colnames `regressors' _cons
	matrix colnames `b_out' = `colnames'
	matrix colnames `V_out' = `colnames'
	matrix rownames `V_out' = `colnames'

	ereturn post `b_out' `V_out', esample(`touse') depname(`depvar')

	ereturn scalar N           = `N_final'
	ereturn scalar converged   = `converged_final'
	ereturn scalar iterations  = `iter_final'
	ereturn scalar deviance    = `dev_final'

	ereturn local cmd          "ivppmlhdfejl"
	ereturn local cmdline      "ivppmlhdfejl `0'"
	ereturn local depvar       "`depvar'"
	ereturn local exogvars     "`exog'"
	ereturn local endogvars    "`endog'"
	ereturn local instruments  "`instruments'"
	ereturn local absorb       "`absorb_raw'"
	ereturn local vcetype      "`vcetype'"
	ereturn local clustvar     "`clustvar'"
	ereturn local title        "IV-PPML with HDFE (Julia)"

	// ---------- 11.  Display ----------
	di as txt _n "{hline 78}"
	di as txt "IV-PPML with High-Dimensional Fixed Effects  (Julia backend v0.2)"
	di as txt "{hline 78}"
	di as txt "Dependent variable: " as res "`depvar'"
	di as txt "Endogenous:         " as res "`endog'"
	di as txt "Instruments:        " as res "`instruments'"
	di as txt "Absorbed FE:        " as res "`absorb_raw'"
	di as txt "Number of obs:      " as res e(N)
	if e(converged) == 1 {
		di as txt "Converged:          " as res "yes" ///
			as txt "  (iterations = " as res e(iterations) as txt ")"
	}
	else {
		di as txt "Converged:          " as err "no" ///
			as txt "  (iterations = " as res e(iterations) as txt ")"
	}
	di as txt "{hline 78}"
	_coef_table, `eform'
	di as txt "{hline 78}"
	di as txt "Endogenous: " as res "`endog'"
	di as txt "Instruments:" as res " `instruments'"
end
