*! ivppmlhdfejl 0.9.4  14apr2026
*! (IV-)PPML with High-Dimensional Fixed Effects — Julia backend
*! Architecture: Path 3 — single Julia call with solver reuse (update_weights!)
*! Mirrors reghdfejl / GLFixedEffectModels.jl pattern
*!
*! Algorithm: IRLS-IV  (Mullahy 1997 GMM via iteratively reweighted 2SLS)
*!            Without IV specification, reduces to standard PPML (ppmlhdfejl)
*! FE absorption via FixedEffects.jl (iterative demeaning)
*! Authors: Ohyun Kwon
*!
*! Syntax (IV):
*!   ivppmlhdfejl depvar [exogvars] (endogvars = excluded_instruments) [if] [in] [pw], ///
*!       absorb(absvars) [options]
*! Syntax (PPML, no IV):
*!   ivppmlhdfejl depvar [regressors] [if] [in] [pw], absorb(absvars) [options]
*!
*! Options: vce(cluster c1 [c2]) tolerance(#) itolerance(#) maxiterations(#) verbose(#)
*!          exposure(varname) offset(varname) d(name) separation(string)
*!          standardize guess(simple|mean)

// Mata helper used by Patch B below: when Julia's v0.9.4 two-stage
// collinearity pipeline drops some columns, this pads b / V back to the
// full expected width, placing zeros at the dropped positions so that
// _ms_findomitted flags them as `o.var` in the coefficient table.
cap mata: mata drop _ivppml_pad_bV()
mata:
void _ivppml_pad_bV(string scalar bname, string scalar Vname,
                    string scalar b_out, string scalar V_out,
                    string rowvector keep)
{
    real rowvector b_got, b_pad
    real matrix    V_got, V_pad
    real scalar    K_exp, K_got, i, j, pi, pj

    b_got = st_matrix(bname)
    V_got = st_matrix(Vname)
    K_exp = cols(keep)
    K_got = cols(b_got)

    b_pad = J(1, K_exp, 0)
    V_pad = J(K_exp, K_exp, 0)

    pi = 0
    for (i = 1; i <= K_exp; i++) {
        if (strtoreal(keep[i]) == 1) {
            pi++
            b_pad[1, i] = b_got[1, pi]
            pj = 0
            for (j = 1; j <= K_exp; j++) {
                if (strtoreal(keep[j]) == 1) {
                    pj++
                    V_pad[i, j] = V_got[pi, pj]
                }
            }
        }
    }
    st_matrix(b_out, b_pad)
    st_matrix(V_out, V_pad)
}
end

program ivppmlhdfejl, eclass
	version 15.0

	if replay() {
		if "`e(cmd)'" != "ivppmlhdfejl" error 301
		// Use a dedicated Replay subroutine so display options (eform, irr,
		// level, noheader, notable) parse correctly instead of being passed
		// raw to _coef_table (which only knows eform("label"), not eform).
		Replay `0'
		exit
	}

	// Cold-start handling: the very first jl call in a new Stata session can
	// fail with a transient initialization error. Make sure ivppmlhdfejl_load
	// runs and the package is loaded BEFORE we try to GetEnv. If GetEnv still
	// fails on the first attempt, retry once.
	cap ivppmlhdfejl_load
	if _rc {
		// Allow the initial install / precompile to surface its real error
		ivppmlhdfejl_load
	}

	cap qui jl GetEnv
	if _rc {
		// One more attempt — first-call init transient
		qui jl GetEnv
	}
	local env `r(env)'

	cap noi _ivppmlhdfejl `0'
	local rc = _rc

	* Restore previous Julia environment
	cap qui jl SetEnv `env'

	exit `rc'
end


// ==========================================================================
// Replay subroutine — parses display options before calling _coef_table
// ==========================================================================

program Replay
	syntax [, EForm IRr noHEADer noTABle Level(cilevel) *]
	if "`irr'" != "" local eform "eform"
	_get_diopts diopts, `options' level(`level')

	if "`header'" == "" {
		di as txt _n "{hline 78}"
		local _is_iv = ("`e(endogvars)'" != "")
		if `_is_iv' {
			di as txt "IV-PPML with High-Dimensional Fixed Effects  (Julia backend)"
		}
		else {
			di as txt "PPML with High-Dimensional Fixed Effects  (Julia backend)"
		}
		di as txt "{hline 78}"
		local K = e(df_m)
		di as txt "Dependent variable: " _col(24) as res "`e(depvar)'" ///
			_col(51) as txt "No. of obs" _col(68) "=" _col(70) as res %10.0gc e(N)
		if `_is_iv' {
			di as txt "Endogenous:         " _col(24) as res "`e(endogvars)'" ///
				_col(51) as txt "Residual df" _col(68) "=" _col(70) as res %10.0gc e(df)
			di as txt "Instruments:        " _col(24) as res "`e(instruments)'" ///
				_col(51) as txt "Wald chi2(`K')" _col(68) "=" _col(70) as res %10.2f e(chi2)
		}
		else {
			di as txt _col(51) "Residual df" _col(68) "=" _col(70) as res %10.0gc e(df)
			di as txt _col(51) "Wald chi2(`K')" _col(68) "=" _col(70) as res %10.2f e(chi2)
		}
		di as txt "Absorbed FE:        " _col(24) as res "`e(absvars)'" ///
			_col(51) as txt "Prob > chi2" _col(68) "=" _col(70) as res %10.4f chi2tail(`K', e(chi2))
		if "`e(offset)'" != "" {
			di as txt "Offset:             " _col(24) as res "`e(offset)'" ///
				_col(51) as txt "Pseudo R2" _col(68) "=" _col(70) as res %10.4f e(r2_p)
		}
		else {
			di as txt _col(51) "Pseudo R2" _col(68) "=" _col(70) as res %10.4f e(r2_p)
		}
		di as txt "Log pseudolikelihood = " as res %12.0g e(ll) ///
			_col(51) as txt "Deviance" _col(68) "=" _col(70) as res %10.4g e(deviance)
		if "`e(vce)'" == "cluster" {
			local ncv = e(N_clustervars)
			forvalues j = 1/`ncv' {
				di _col(51) as txt "No. of clusters" _col(68) "=" _col(70) as res %10.0fc e(N_clust`j') ///
					as txt "  (`e(clustvar`j')')"
			}
		}
		if "`e(converged)'" == "1" {
			di as txt "Converged:          " as res "yes" ///
				as txt "  (iterations = " as res e(iterations) as txt ")"
		}
		else {
			di as txt "Converged:          " as err "no" ///
				as txt "  (iterations = " as res e(iterations) as txt ")"
		}
		if e(num_singletons) > 0 {
			di as txt "(`=e(num_singletons)' singleton observations dropped)"
		}
		if e(num_separated) > 0 {
			di as txt "(`=e(num_separated)' separated observations dropped)"
		}
	}
	if "`table'" == "" {
		di as txt "{hline 78}"
		_coef_table, `eform' `diopts'
		di as txt "{hline 78}"
		if "`e(endogvars)'" != "" {
			di as txt "Endogenous: " as res "`e(endogvars)'"
			di as txt "Instruments:" as res " `e(instruments)'"
		}
	}
end


// ==========================================================================
// Parse syntax:  depvar [exog] [(endog = instruments)]  [if] [in], options
// Without parentheses: all regressors are exogenous (standard PPML)
// ==========================================================================

program ParseIV, sclass
	sreturn clear
	local cmd `"`0'"'

	// Detect IV vs PPML: check if ( appears before the first comma
	// In IV syntax: depvar [exog] (endog = inst) [if] [in], options
	// In PPML syntax: depvar [regressors] [if] [in], options
	local cpos   = strpos(`"`cmd'"', ",")
	local popen  = strpos(`"`cmd'"', "(")

	// IV mode if ( exists AND comes before the comma (or no comma)
	local is_iv = 0
	if `popen' > 0 {
		if `cpos' == 0 | `popen' < `cpos' {
			local pclose = strpos(`"`cmd'"', ")")
			if `pclose' > `popen' {
				local is_iv = 1
			}
		}
	}

	if !`is_iv' {
		// No IV specification — standard PPML
		// Syntax: depvar [regressors] [if] [in] , options
		// Strategy: extract just depvar, then pass everything else
		// (including if/in/options) to _ivppmlhdfejl's syntax parser
		gettoken depvar rest : cmd
		// rest = "x1 x2 [if ...] [in ...], absorb(...) ..."
		// We need to split exog vars from [if] [in] , options
		// Find the comma — everything before it that isn't if/in is exog
		local rest = strtrim(`"`rest'"')
		local cpos2 = strpos(`"`rest'"', ",")
		if `cpos2' == 0 {
			di as err "absorb() is required"
			exit 198
		}
		local vpart = strtrim(substr(`"`rest'"', 1, `cpos2' - 1))
		local after = substr(`"`rest'"', `cpos2', .)

		// vpart may contain: "x1 x2" or "x1 x2 if x1>0" or "x1 x2 in 1/100"
		// or "x1 x2 [pw=wt]"; parse out [if] [in] [weight].
		// `fv ts` lets users pass factor / time-series operators like
		// `i.region`, `c.x##i.year`, `L.x` in PPML mode (IV mode already
		// works via string splicing).
		local 0 `"`vpart'"'
		syntax [varlist(fv ts default=none)] [if] [in] [pw fw/]
		local exog `varlist'
		local _wexp ""
		if "`weight'" != "" local _wexp `"[`weight'=`exp']"'
		// Reconstruct after with if/in/weight prepended
		local after `"`if' `in' `_wexp' `after'"'

		sreturn local depvar       `"`depvar'"'
		sreturn local exog         `"`exog'"'
		sreturn local endog        ""
		sreturn local instruments  ""
		sreturn local after        `"`after'"'
		exit
	}

	// IV specification present
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

	// Save full command line before `0' gets rewritten by ParseIV / syntax.
	local full_cmdline `"`0'"'

	// ---------- 1.  Parse IV specification ----------
	ParseIV `0'
	local depvar       `s(depvar)'
	local exog         `s(exog)'
	local endog        `s(endog)'
	local instruments  `s(instruments)'
	local after_paren  `s(after)'

	// ---------- 2.  Parse options ----------
	local 0 `"`after_paren'"'
	syntax [if] [in] [pw fw/] , [Absorb(string) NOAbsorb ///
		VCE(string) CLuster(string) ///
		 TOLerance(real 1e-8) ITOLerance(real -1) MAXITerations(integer 1000) ///
		 Verbose(integer 0) noLOG EForm IRr ///
		 EXPosure(varname) OFFset(varname) ///
		 D(name) SEParation(string) ///
		 STANDardize GUESS(string) ///
		 gpu THReads(integer 0) KEEPSINgletons compact]

	// absorb / noabsorb (mutually exclusive; one is required)
	if "`absorb'" == "" & "`noabsorb'" == "" {
		di as err "must specify either absorb() or noabsorb"
		exit 198
	}
	if "`absorb'" != "" & "`noabsorb'" != "" {
		di as err "cannot specify both absorb() and noabsorb"
		exit 198
	}
	if "`noabsorb'" != "" {
		tempvar _one
		gen byte `_one' = 1
		local absorb "`_one'"
	}

	// Detect slope-only absorb (e.g., `absorb(id#c.t)`). The Julia backend
	// silently drops the obs count to 0 on such specs because the FE has
	// no intercept to absorb. Refuse cleanly; `##c.` expands to both the
	// intercept and slope FE and is safe. Mirrors native ivppmlhdfe.
	foreach _abs_tok of local absorb {
		if regexm("`_abs_tok'", "#c\.") & !regexm("`_abs_tok'", "##c\.") {
			di as err "slope-only absorb (`_abs_tok') not supported by ivppmlhdfejl; use `##c.` to include the intercept FE"
			exit 198
		}
	}

	// Validate guess()
	if !inlist("`guess'", "", "simple", "mean", "default") {
		di as err "guess() must be simple, mean, or default"
		exit 198
	}
	local guess_jl ":default"
	if "`guess'" == "mean"   local guess_jl ":mean"
	if "`guess'" == "simple" local guess_jl ":default"

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

	// d() option — auto-name to _ivppmlhdfejl_d so `predict mu` works without
	// the user having to remember to ask for d() at estimation. Drop any
	// existing copy so re-fits in loops don't error rc=110. Mirrors native.
	if "`d'" == "" {
		local d "_ivppmlhdfejl_d"
	}
	capture drop `d'

	// separation (canonicalise + validate against allowed list)
	if inlist("`separation'", "", "def", "default", "on", "auto", "standard") {
		local separation "fe"
	}
	if inlist("`separation'", "all", "full") {
		local separation "fe"
		di as txt "(note: ivppmlhdfejl bridge only implements separation(fe); simplex/relu/mu are native-only)"
	}
	if inlist("`separation'", "no", "off", "none") {
		local separation ""
	}
	// Anything other than fe or empty is rejected — bridge only does fe
	local _allowed_sep "fe"
	if "`separation'" != "" & "`: list separation - _allowed_sep'" != "" {
		di as err "ivppmlhdfejl: separation(`separation') not supported in the Julia bridge"
		di as err "  valid values: default | fe | none"
		di as err "  (use the native ivppmlhdfe command for simplex / relu / mu)"
		exit 198
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
	// Use fvrevar to materialise factor / time-series operators (e.g.,
	// i.exp -> set of byte tempvars). unab alone cannot expand i.x.
	// fvrevar may emit base-level dummies that are constant — those would
	// trip the Julia engine's zero-variance guard, so filter them locally.
	// Helper logic: fvrevar emits base-level dummies as constant tempvars
	// (e.g., the omitted base of `i.year`). Filter those out so the Julia
	// engine's zero-variance guard doesn't error. Apply to exog, endog, and
	// instruments equally.
	//
	// To get FACTOR-AWARE display names for e(b) (so the coef table shows
	// `2.region` instead of `__000001`), we call fvrevar twice: once with
	// `, list` to get the display labels, once normally to materialise the
	// tempvars. The two lists are the SAME length and in the SAME order, so
	// we can zip them into `exog` / `_exog_labels`.
	// fvexpand gives display labels with `Nb.var` base markers, which is
	// parallel (same order, same length) to the tempvars fvrevar produces.
	// Zip them and drop: (a) base-level dummies (`Nb.var`), (b) zero-variance
	// columns (shouldn't happen after fvexpand but defensive).
	local _exog_labels ""
	local _endog_labels ""
	local _inst_labels  ""
	foreach _role in exog endog instruments {
		if "``_role''" == "" continue
		fvexpand ``_role''
		local _lbl_list `r(varlist)'
		fvrevar ``_role''
		local _tvar_list `r(varlist)'
		local _kept ""
		local _klbl ""
		local _zerovar ""
		local i_lbl 0
		foreach _t of local _tvar_list {
			local ++i_lbl
			local _lbl : word `i_lbl' of `_lbl_list'
			if regexm("`_lbl'", "^[0-9]+b[n]?\.") continue
			capture confirm variable `_t'
			if _rc continue
			quietly summarize `_t', meanonly
			if r(min) == r(max) {
				local _zerovar `_zerovar' `_lbl'
				continue
			}
			local _kept `_kept' `_t'
			local _klbl `_klbl' `_lbl'
		}
		// Surface zero-variance instruments / regressors with a clear error
		// instead of silently dropping them and reporting "not identified" or
		// missing-coefficient downstream. Mirrors native ivppmlhdfe.
		if "`_zerovar'" != "" {
			if "`_role'" == "instruments" {
				di as err "instrument(s) with zero variance (constant column): `_zerovar'"
				exit 459
			}
			else {
				di as err "regressor(s) with zero variance (constant column): `_zerovar'"
				exit 459
			}
		}
		local `_role' `_kept'
		if "`_role'" == "exog"        local _exog_labels  `_klbl'
		if "`_role'" == "endog"       local _endog_labels `_klbl'
		if "`_role'" == "instruments" local _inst_labels  `_klbl'
	}
	if "`clustvars'" != "" {
		local clustvars_unab ""
		local clustvars_display ""
		foreach cv of local clustvars {
			// Handle interaction notation (e.g., c1#c2 -> egen group)
			if strpos("`cv'", "#") > 0 {
				local cv_clean : subinstr local cv "##" "#", all
				local cv_clean : subinstr local cv_clean "#" " ", all
				local cv_clean : subinstr local cv_clean "i." "", all
				local cv_clean : subinstr local cv_clean "c." "", all
				tempvar _bridge_clust_int
				quietly egen long `_bridge_clust_int' = group(`cv_clean')
				local clustvars_unab `clustvars_unab' `_bridge_clust_int'
				local clustvars_display `clustvars_display' `cv'
			}
			else {
				// Handle string clusters by encoding to numeric (matches native).
				unab cv : `cv'
				capture confirm string variable `cv'
				if !_rc {
					tempvar _bridge_clust_str
					quietly egen long `_bridge_clust_str' = group(`cv')
					local clustvars_unab `clustvars_unab' `_bridge_clust_str'
					local clustvars_display `clustvars_display' `cv'
				}
				else {
					local clustvars_unab `clustvars_unab' `cv'
					local clustvars_display `clustvars_display' `cv'
				}
			}
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

	// Validate exposure > 0 BEFORE markout drops missing-offset obs.
	// `offset' is the tempvar `offset_var' = ln(`exposure'), missing when
	// exposure <= 0; without this preflight markout would silently drop
	// the bad rows and the user would see an N off by a few. Mirrors native.
	if "`exposure'" != "" {
		quietly count if `exposure' <= 0 & `touse'
		if r(N) > 0 {
			di as err "exposure() must be greater than zero"
			exit 459
		}
	}
	// Reject negative weights (mirrors native).
	if "`weight'" != "" {
		quietly count if `exp' < 0 & `touse'
		if r(N) > 0 {
			di as err "negative weights encountered"
			exit 402
		}
		// Drop zero-weight rows — markout only drops missing. Zero-weight rows
		// contribute nothing and can crash the Julia `st_store`/partial-out
		// path with a Mata conformability error when they flow through.
		quietly replace `touse' = 0 if `exp' == 0 & `touse'
	}

	if "`offset'" != "" markout `touse' `offset'

	quietly count if `depvar' < 0 & `touse'
	if r(N) > 0 {
		di as err "`depvar' must be nonnegative"
		exit 459
	}

	quietly count if `touse'
	local N_full = r(N)
	if `N_full' == 0 error 2000

	// Pre-flight degenerate-input check: bail before crossing into Julia,
	// which would otherwise leak a SingularException / FieldError stack
	// trace on tiny / under-identified samples.
	local _Kpre = `n_exog' + `n_endog' + 1
	if `N_full' <= `_Kpre' {
		di as err "ivppmlhdfejl: insufficient observations (N=`N_full', K+1=`_Kpre')"
		exit 2001
	}

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

	// ---------- 4b.  Singleton drop (default; suppressed by keepsingletons) ----------
	// Iteratively drops obs that are alone in any FE group, matching the
	// reghdfe / native ivppmlhdfe convention. Without this, the bridge's
	// e(N), e(df), and _cons drift away from the native command.
	// Build the cleaned varlist for each absorb token (JOINT interaction
	// spec) so singletons are counted on the interaction (exp#year) not on
	// individual vars. We build the raw var-lists ONCE, then rebuild the
	// group ids INSIDE the cascade loop so counts are fresh on each pass.
	local _sing_spec_list ""
	local _sing_spec_count 0
	foreach avar of local absorb {
		local avar_s "`avar'"
		local avar_s : subinstr local avar_s "##" "#", all
		local avar_s : subinstr local avar_s "#" " ", all
		local avar_s : subinstr local avar_s "i." "", all
		local avar_s : subinstr local avar_s "c." "", all
		local _ok 1
		foreach v of local avar_s {
			capture confirm numeric variable `v'
			if _rc {
				local _ok 0
			}
		}
		if `_ok' {
			local ++_sing_spec_count
			local _sing_spec_`_sing_spec_count' "`avar_s'"
		}
	}

	local n_singletons = 0
	if "`keepsingletons'" == "" {
		local sing_changed = 1
		while `sing_changed' {
			local sing_changed = 0
			forvalues _ii = 1/`_sing_spec_count' {
				tempvar _sgrp grpcount
				// Rebuild the joint group id on the CURRENT (shrunken) touse.
				// Without this, after the first pass drops obs from one token,
				// the stale group id for another token still shows groups of
				// size >1 and the cascade fails to converge.
				quietly egen long `_sgrp' = group(`_sing_spec_`_ii'') if `touse'
				quietly bysort `_sgrp': gen long `grpcount' = _N if `touse'
				quietly count if `grpcount' == 1 & `touse'
				if r(N) > 0 {
					quietly replace `touse' = 0 if `grpcount' == 1 & `touse'
					local n_singletons = `n_singletons' + r(N)
					local sing_changed = 1
				}
				drop `_sgrp' `grpcount'
			}
		}
		quietly count if `touse'
		if r(N) == 0 {
			di as err "no observations remaining after singleton drop (try keepsingletons)"
			exit 2000
		}
	}

	// ---------- 5.  d() variable ----------
	if "`d'" != "" {
		quietly gen double `d' = .
	}

	// ---------- 6.  Parse absorb ----------
	// Each token can be a single var (id), an i.var, or an interaction
	// using # / ## (e.g., exp#year, i.exp#i.year, c.x##i.year). For each
	// token we strip i. / c. operators, expand # to spaces, and feed the
	// resulting var list to egen group() to build a numeric FE id. This
	// matches what reghdfe does internally for these patterns.
	local absorb_raw `"`absorb'"'

	local fe_vars ""
	local fe_syms ""
	local n_fe 0
	foreach avar of local absorb {
		local avar_clean "`avar'"
		// "##" produces both main effects and interaction. We approximate
		// it as a pure interaction (the main effects would otherwise be
		// absorbed by the interaction anyway).
		local avar_clean : subinstr local avar_clean "##" "#", all
		// Replace # with space so egen group() sees a varlist
		local avar_clean : subinstr local avar_clean "#" " ", all
		// Strip factor / time-series operators
		local avar_clean : subinstr local avar_clean "i." "", all
		local avar_clean : subinstr local avar_clean "c." "", all
		local ++n_fe
		tempvar fe_`n_fe'
		capture qui egen long `fe_`n_fe'' = group(`avar_clean') if `touse'
		if _rc {
			di as err "ivppmlhdfejl: cannot parse absorb token `avar'"
			di as err "  (cleaned to: `avar_clean')"
			exit 198
		}
		local fe_vars `fe_vars' `fe_`n_fe''
		local fe_syms `fe_syms' `fe_`n_fe''
	}

	// Compute an approximate df_a (absorbed degrees of freedom) from the
	// FE tempvars, so the residual df printed in the coef table isn't
	// wildly off for unclustered fits. Formula: sum(distinct levels) -
	// n_fe + 1 (one FE keeps its constant). Matches reghdfe's default
	// heuristic in the common non-nested case.
	local df_a = 0
	forvalues _j = 1/`n_fe' {
		quietly tabulate `fe_`_j'' if `touse'
		local df_a = `df_a' + r(r)
	}
	if `n_fe' > 0 {
		local df_a = `df_a' - `n_fe' + 1
	}
	if `df_a' < 0 local df_a = 0

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
	if `n_endog' > 0 {
		foreach v of local endog {
			if "`jl_endog'" == "" local jl_endog Symbol("`v'")
			else                 local jl_endog `jl_endog', Symbol("`v'")
		}
	}
	local jl_inst ""
	if `n_inst' > 0 {
		foreach v of local instruments {
			if "`jl_inst'" == "" local jl_inst Symbol("`v'")
			else                local jl_inst `jl_inst', Symbol("`v'")
		}
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

	// Build new-options keyword fragment
	local jl_stdopt ""
	if "`standardize'" != "" local jl_stdopt , standardize = true
	local jl_guessopt , guess = `guess_jl'
	local jl_itolopt ""
	if `itolerance' > 0 local jl_itolopt , itol = `itolerance'

	// Estimate!
	// v0.9.4: wrap the Julia call in try/catch so runaway-divergence /
	// rank-deficient / non-convergence errors from ivppml_reg surface as
	// clean Stata rc codes (9003 / 430 / 9004) instead of leaking a Julia
	// stacktrace and then throwing rc=111 on the first post-estimation
	// tempvar access. Communication via st_numscalar (st_local is not
	// available in all jl package versions).
	tempname t_jl_rc
	_jl: global ivppml_jl_rc = 0.0;
	_jl: try; global m = IVPPMLFixedEffectModels.ivppml_reg(df, ///
		Symbol("`depvar'"), ///
		Symbol[`jl_exog'], ///
		Symbol[`jl_endog'], ///
		Symbol[`jl_inst'], ///
		Symbol[`jl_fes'], ///
		tol = `tolerance', maxiter = `maxiterations', verbose = `verbose' ///
		`jl_clustopt' `jl_wtopt' `jl_offopt' ///
		`jl_stdopt' `jl_guessopt' `jl_itolopt'); ///
		catch _e; ///
		global m = nothing; ///
		_msg = sprint(showerror, _e); ///
		if occursin("rc=9003", _msg); global ivppml_jl_rc = 9003.0; ///
		elseif occursin("rc=430", _msg); global ivppml_jl_rc = 430.0; ///
		elseif occursin("2SLS solve failed", _msg) || occursin("rank-deficient", _msg); global ivppml_jl_rc = 9004.0; ///
		else; global ivppml_jl_rc = 9000.0; end; ///
		println("ivppmlhdfejl: ", _msg); ///
		end;

	// Free DataFrame regardless of success/failure
	_jl: df = nothing;

	// Check if Julia errored
	_jl: st_numscalar("`t_jl_rc'", ivppml_jl_rc);
	local jl_rc = `t_jl_rc'
	if `jl_rc' != 0 {
		di as err "ivppmlhdfejl: Julia estimation failed (rc=`jl_rc')"
		di as err "  See the Julia diagnostic message above."
		_jl: global ivppml_jl_rc = 0.0; global m = nothing;
		exit `jl_rc'
	}

	// ---------- 12.  Retrieve results ----------
	tempname b_out V_out t

	_jl: st_numscalar("`t'", m.N);
	local N_final = `t'

	_jl: st_numscalar("`t'", m.converged ? 1.0 : 0.0);
	local converged_final = `t'

	// Defence-in-depth: Julia's v0.9.3 guards now hard-error on divergence
	// and non-convergence, so this branch is rarely hit. Kept in case a
	// future code path slips `converged=false` through without raising.
	if `converged_final' == 0 {
		di as err "ivppmlhdfejl: Julia engine did not converge"
		di as err "  coefficients are numerically meaningless and will not be posted."
		di as err "  Try standardize, mu_separation, or a simpler FE structure."
		_jl: m = nothing; global ivppml_m = nothing;
		exit 9003
	}

	_jl: st_numscalar("`t'", Float64(m.iterations));
	local iter_final = `t'

	_jl: st_numscalar("`t'", m.deviance);
	local dev_final = `t'

	_jl: st_numscalar("`t'", m.ll);
	local ll_final = `t'

	_jl: st_numscalar("`t'", m.ll_0);
	local ll0_final = `t'

	_jl: st_numscalar("`t'", Float64(m.rank));
	local rank_final = `t'

	_jl: st_numscalar("`t'", Float64(m.num_sep_mu));
	local nsep_mu_final = `t'

	// Flag rank deficiency so users see a clear warning rather than wonder
	// why one factor-variable dummy has a huge SE. The Julia engine reports
	// rank < K when any column of Xhat'WX is collinear (typically a factor
	// dummy nested in an absorbed FE dimension).
	local K_total = `n_exog' + `n_endog'
	if `rank_final' < `K_total' {
		di as err "(warning: design matrix rank `rank_final' < K=`K_total'; " ///
			"some regressors are collinear with the absorbed FE.)"
		di as err "  Suspect: factor-variable dummies nested in absorb(). " ///
			"The Julia engine keeps all columns in e(b) but their coefficients are not identified."
	}

	// Get coefficient vector
	_jl: _ivppml_b = collect(IVPPMLFixedEffectModels.coef(m)');
	jl GetMatFromMat `b_out', source(_ivppml_b)

	// Get VCE matrix
	_jl: _ivppml_V = IVPPMLFixedEffectModels.vcov(m);
	jl GetMatFromMat `V_out', source(_ivppml_V)

	// Get per-block drop counts from Julia (st_numscalar always works).
	tempname t_ndrop_exog t_ndrop_endog t_ndrop_inst
	_jl: st_numscalar("`t_ndrop_exog'", Float64(m.n_dropped_exog));
	_jl: st_numscalar("`t_ndrop_endog'", Float64(m.n_dropped_endog));
	_jl: st_numscalar("`t_ndrop_inst'", Float64(m.n_dropped_inst));
	local ndrop_exog  = `t_ndrop_exog'
	local ndrop_endog = `t_ndrop_endog'
	local ndrop_inst  = `t_ndrop_inst'

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

	// ---------- 13.  Post results ----------
	// Use the factor-aware display labels captured from `fvrevar , list`
	// so e(b) shows `2.region` etc. instead of the raw `__000001` tempvar.
	local colnames `_exog_labels' `_endog_labels' _cons
	local K_expected : word count `colnames'
	local K_got = colsof(`b_out')

	// v0.9.4: Julia's two-stage collinearity pipeline may have dropped some
	// columns (fully absorbed by FE, or collinear with other regressors /
	// instruments). In that case `b_out` is shorter than the expected Stata
	// colnames. Pad with zeros at the dropped positions so e(b) / e(V)
	// match the Stata-native "omitted" display convention.
	if `K_got' < `K_expected' {
		local n_drop = `K_expected' - `K_got'
		di as txt "(two-stage collinearity: Julia dropped `n_drop' column(s); " ///
			"they appear as omitted in e(b).)"

		// Build keepvec using per-block drop counts from Julia.
		// Layout: [exog_1..exog_n_exog] [endog_1..endog_n_endog] [_cons]
		// For the common case (all exog absorbed by FE, no endog dropped),
		// ndrop_exog == n_exog and ndrop_endog == 0.
		local keepvec ""
		local n_exog_kept = `n_exog' - `ndrop_exog'
		local n_endog_kept = `n_endog' - `ndrop_endog'
		// Exog: first n_exog_kept are kept, rest dropped
		forvalues j = 1/`n_exog' {
			if `j' <= `n_exog_kept' {
				local keepvec `keepvec' 1
			}
			else {
				local keepvec `keepvec' 0
			}
		}
		// Endog: first n_endog_kept are kept, rest dropped
		forvalues j = 1/`n_endog' {
			if `j' <= `n_endog_kept' {
				local keepvec `keepvec' 1
			}
			else {
				local keepvec `keepvec' 0
			}
		}
		// _cons always kept
		local keepvec `keepvec' 1

		tempname b_pad V_pad
		mata: _ivppml_pad_bV("`b_out'", "`V_out'", "`b_pad'", "`V_pad'", ///
			tokens("`keepvec'"))
		matrix `b_out' = `b_pad'
		matrix `V_out' = `V_pad'
	}

	// Clean up Julia memory
	_jl: _ivppml_b = _ivppml_V = nothing;
	_jl: m = nothing; global ivppml_m = nothing;

	matrix colnames `b_out' = `colnames'
	matrix colnames `V_out' = `colnames'
	matrix rownames `V_out' = `colnames'

	// Mark omitted (collinear) factor-variable dummies before posting so the
	// coefficient table shows them as `o.` instead of `__000001` with absurd
	// SEs. Native ivppmlhdfe does the same call after assembling b/V.
	// _ms_findomitted flags columns where b==0 AND diag(V)==0 — exactly
	// the shape the padding branch above produces for dropped columns.
	capture _ms_findomitted `b_out' `V_out'

	ereturn post `b_out' `V_out', esample(`touse') depname(`depvar')

	// Scalars
	ereturn scalar N           = `N_final'
	ereturn scalar converged   = `converged_final'
	ereturn scalar iterations  = `iter_final'
	// Alias e(ic) — Stata-standard iteration counter; many downstream tools
	// (estat ic, esttab) read e(ic) rather than e(iterations).
	ereturn scalar ic          = `iter_final'
	ereturn scalar deviance    = `dev_final'
	ereturn scalar ll          = `ll_final'
	ereturn scalar ll_0        = `ll0_final'
	ereturn scalar r2_p        = 1 - `ll_final' / `ll0_final'
	ereturn scalar N_full      = `N_full'
	ereturn scalar num_separated = `n_separated'
	ereturn scalar num_singletons = `n_singletons'

	local K = `n_exog' + `n_endog'
	ereturn scalar df_m        = `rank_final'
	ereturn scalar rank        = `rank_final'
	ereturn scalar num_sep_mu  = `nsep_mu_final'

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
		// Warn on a degenerate G=1 cluster — variance is not identified.
		// Native ivppmlhdfe falls back to robust automatically; here we
		// match by warning loudly so users do not silently trust a label.
		forvalues j = 1/`n_clustvars' {
			if `G`j'' <= 1 {
				di as err "warning: cluster variable #`j' has only `G`j'' cluster(s); cluster VCE not defined"
				di as err "         Stata's _coef_table will still print 'cluster' but the SEs are unreliable."
			}
		}
		ereturn scalar N_clustervars = `n_clustvars'
		// Use the user-facing display names (clustvars_display) if we built
		// any tempvars for string / interaction clusters above; otherwise
		// fall back to clustvars.
		local _show "`clustvars_display'"
		if "`_show'" == "" local _show "`clustvars'"
		forvalues j = 1/`n_clustvars' {
			local cvar : word `j' of `_show'
			ereturn local clustvar`j' "`cvar'"
			ereturn scalar N_clust`j' = `G`j''
		}
		// e(df) = min clusters - 1
		local min_G = `G1'
		if `n_clustvars' > 1 {
			forvalues j = 2/`n_clustvars' {
				if `G`j'' < `min_G' local min_G = `G`j''
			}
		}
		ereturn scalar df      = `min_G' - 1
		ereturn scalar N_clust = `min_G'
		// e(clustvar) is the joint name string read by _coef_table when it
		// prints the "Std. err. adjusted for X clusters in ..." header.
		ereturn local clustvar "`_show'"
	}
	else {
		// Unclustered: residual df accounts for the absorbed FE coefficients.
		ereturn scalar df   = `N_final' - `K' - `df_a'
		ereturn scalar df_a = `df_a'
	}

	// Locals
	local is_iv = (`n_endog' > 0)

	ereturn local cmd          "ivppmlhdfejl"
	ereturn local cmdline      `"ivppmlhdfejl `full_cmdline'"'
	ereturn local depvar       "`depvar'"
	ereturn local exogvars     "`exog'"
	ereturn local indepvars    "`regressors' _cons"
	if "`weight'" != "" {
		ereturn local wtype    "`weight'"
		ereturn local wexp     "= `exp'"
	}
	if `is_iv' {
		ereturn local endogvars    "`endog'"
		ereturn local instruments  "`instruments'"
	}
	ereturn local absorb       "`absorb_raw'"
	ereturn local absvars      "`absorb_raw'"
	// Stata convention: e(vce) is the lowercase internal name, e(vcetype)
	// is the display label used by _coef_table ("Robust" with capital R).
	ereturn local vce          "`vcetype'"
	ereturn local vcetype      "Robust"
	// Use display names for e(clustvar), not the encoded tempvar list — _coef_table
	// reads e(clustvar) for the "Std. err. adjusted for ..." annotation.
	if "`clustvars_display'" != "" {
		ereturn local clustvar "`clustvars_display'"
	}
	else {
		ereturn local clustvar "`clustvars'"
	}
	if `is_iv' {
		ereturn local title    "IV-PPML with HDFE (Julia)"
	}
	else {
		ereturn local title    "PPML with HDFE (Julia)"
	}
	ereturn local chi2type     "Wald"
	ereturn local predict      "ivppmlhdfe_p"
	ereturn local marginsok    "default"
	ereturn local properties   "b V"
	ereturn local separation   "`separation'"
	if "`offset_display'" != "" ereturn local offset "`offset_display'"
	if "`d'" != "" ereturn local d "`d'"

	// Populate e(rss) / e(rmse) for esttab / estat ic compatibility (mirror
	// native ivppmlhdfe). Must run AFTER `e(predict) = "ivppmlhdfe_p"` so
	// `predict mu` dispatches to our handler.
	tempvar _muh _resh
	capture noisily predict double `_muh' if e(sample), mu
	if !_rc {
		quietly gen double `_resh' = (`depvar' - `_muh')^2 if e(sample)
		quietly summarize `_resh' if e(sample), meanonly
		if r(N) > 0 {
			ereturn scalar rss = r(sum)
			if e(df) > 0 ereturn scalar rmse = sqrt(r(sum) / e(df))
		}
	}
	capture drop `_muh' `_resh'

	// ---------- 14.  Display ----------
	di as txt _n "{hline 78}"
	if `is_iv' {
		di as txt "IV-PPML with High-Dimensional Fixed Effects  (Julia backend v0.9.3)"
	}
	else {
		di as txt "PPML with High-Dimensional Fixed Effects  (Julia backend v0.9.3)"
	}
	di as txt "{hline 78}"

	* Left-right header (following ppmlhdfe_header / ivppmlhdfe)
	di as txt "Dependent variable: " _col(24) as res "`depvar'" ///
		_col(51) as txt "No. of obs" _col(68) "=" _col(70) as res %10.0gc e(N)
	if `is_iv' {
		di as txt "Endogenous:         " _col(24) as res "`endog'" ///
			_col(51) as txt "Residual df" _col(68) "=" _col(70) as res %10.0gc e(df)
		di as txt "Instruments:        " _col(24) as res "`instruments'" ///
			_col(51) as txt "Wald chi2(`K')" _col(68) "=" _col(70) as res %10.2f e(chi2)
	}
	else {
		di as txt _col(51) "Residual df" _col(68) "=" _col(70) as res %10.0gc e(df)
		di as txt _col(51) "Wald chi2(`K')" _col(68) "=" _col(70) as res %10.2f e(chi2)
	}
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
			// Use display names (clustvars_display), NOT clustvars which has
			// been rewritten to encoded tempvars like __000000 for string and
			// interaction clusters.
			local cvar : word `j' of `clustvars_display'
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
	if `is_iv' {
		di as txt "Endogenous: " as res "`endog'"
		di as txt "Instruments:" as res " `instruments'"
	}
end
