*! ivppmlhdfe 0.9.4  23apr2026
*! IV-PPML with High-Dimensional Fixed Effects
*! Algorithm: IRLS-IV targeting the additive moment E[q(y-mu)]=0, q=(x',z')',
*!            following Windmeijer & Santos Silva (1997). Solved via iteratively
*!            reweighted 2SLS (Correia, Guimaraes & Zylkin 2020), with FE
*!            concentrated out at each iteration through reghdfe.
*! Authors: Ohyun Kwon, Mario Larch, Jangsu Yoon, Yoto V. Yotov
*! Contact: Ohyun Kwon, theekwonomist@gmail.com
*!
*! Syntax:
*!   ivppmlhdfe depvar [exogvars] (endogvars = excluded_instruments) [if] [in] [pw], ///
*!       absorb(absvars) [vce(cluster c1 [c2]) tolerance(#) maxiterations(#) verbose(#) ///
*!       exposure(varname) offset(varname) d(name) separation(string)]

program ivppmlhdfe, eclass
	version 14.0

	if replay() {
		if "`e(cmd)'" != "ivppmlhdfe" error 301
		Replay `0'
		exit
	}

	capture noisily Estimate `0'
	local rc = c(rc)
	// Cleanup stale Mata externals on error
	if `rc' {
		capture mata: mata drop ivppmlhdfe_HDFE
	}
	exit `rc'
end


program Replay
	syntax [, EForm IRr noHEADer noTABle *]
	if "`irr'" != "" local eform "eform"
	_get_diopts diopts, `options'

	if "`header'" == "" {
		di as txt _n "{hline 78}"
		di as txt "IV-PPML with High-Dimensional Fixed Effects"
		di as txt "{hline 78}"
		local K = e(df_m)
		di as txt "Dependent variable: " _col(24) as res "`e(depvar)'" ///
			_col(51) as txt "No. of obs" _col(68) "=" _col(70) as res %10.0gc e(N)
		di as txt "Endogenous:         " _col(24) as res "`e(endogvars)'" ///
			_col(51) as txt "Residual df" _col(68) "=" _col(70) as res %10.0gc e(df)
		di as txt "Instruments:        " _col(24) as res "`e(instruments)'" ///
			_col(51) as txt "Wald chi2(`K')" _col(68) "=" _col(70) as res %10.2f e(chi2)
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
		if e(converged) == 1 {
			di as txt "Converged:          " as res "yes" ///
				as txt "  (iterations = " as res e(ic) as txt ")"
		}
		else {
			di as txt "Converged:          " as err "no" ///
				as txt "  (iterations = " as res e(ic) as txt ")"
		}
		if e(num_singletons) > 0 {
			di as txt "(`=e(num_singletons)' singleton observations dropped)"
		}
		if e(num_separated) > 0 {
			if e(num_sep_advanced) > 0 | e(num_sep_mu) > 0 {
				local sep_detail ""
				if e(num_sep_fe) > 0 local sep_detail "`=e(num_sep_fe)' by fe"
				if e(num_sep_advanced) > 0 {
					if "`sep_detail'" != "" local sep_detail "`sep_detail', "
					local sep_detail "`sep_detail'`=e(num_sep_advanced)' by simplex/relu"
				}
				if e(num_sep_mu) > 0 {
					if "`sep_detail'" != "" local sep_detail "`sep_detail', "
					local sep_detail "`sep_detail'`=e(num_sep_mu)' by mu"
				}
				di as txt "(`=e(num_separated)' separated observations dropped: `sep_detail')"
			}
			else {
				di as txt "(`=e(num_separated)' separated observations dropped)"
			}
		}
	}

	if "`table'" == "" {
		di as txt "{hline 78}"
		_coef_table, `eform' `diopts'
		di as txt "{hline 78}"
		di as txt "Endogenous: " as res "`e(endogvars)'"
		di as txt "Instruments:" as res " `e(instruments)'"
	}

	// Show absorbed FE footnote (matching ppmlhdfe)
	capture reghdfe_footnote
	// return add  // not needed in eclass context
end


// ==========================================================================
// Parse IV syntax:  depvar [exog] (endog = instruments)  [if] [in], options
// ==========================================================================

program ParseIV, sclass
	sreturn clear
	local cmd `"`0'"'

	// Detect IV vs PPML mode: IV if `(` appears before the first `,`.
	local cpos   = strpos(`"`cmd'"', ",")
	local popen  = strpos(`"`cmd'"', "(")

	local is_iv = 0
	if `popen' > 0 {
		if `cpos' == 0 | `popen' < `cpos' {
			local pclose = strpos(`"`cmd'"', ")")
			if `pclose' > `popen' local is_iv = 1
		}
	}

	if !`is_iv' {
		// Plain PPML mode: depvar [regressors] [if] [in] [pw fw], absorb(...) opts
		// Hand off the entire post-depvar string to the option parser.
		// Use parse(" ,") so a comma immediately after depvar (e.g.,
		// `ivppmlhdfe y, noabsorb` with no regressors) is NOT swallowed
		// into depvar.
		gettoken depvar rest : cmd, parse(" ,")
		local rest = strtrim(`"`rest'"')
		local cpos2 = strpos(`"`rest'"', ",")
		if `cpos2' == 0 {
			di as err "absorb() or noabsorb is required"
			exit 198
		}
		local vpart = strtrim(substr(`"`rest'"', 1, `cpos2' - 1))
		local after = substr(`"`rest'"', `cpos2', .)
		// Split vpart into exog vs [if] [in] [pw fw=...]
		local 0 `"`vpart'"'
		syntax [varlist(default=none fv ts)] [if] [in] [pw fw/]
		local exog `varlist'
		// Re-attach [if] [in] [weight] to after so the main parser sees them.
		local _wexp ""
		if "`weight'" != "" local _wexp `"[`weight'=`exp']"'
		local after `"`if' `in' `_wexp' `after'"'

		sreturn local depvar       `"`depvar'"'
		sreturn local exog         `"`exog'"'
		sreturn local endog        ""
		sreturn local instruments  ""
		sreturn local after        `"`after'"'
		exit
	}

	// IV mode: depvar [exog] (endog = instruments) [if] [in], options
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

	// Save full command line before 0 gets overwritten
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
		TOLerance(real 1e-8) ITOLerance(real -1) ///
		MAXITerations(integer 10000) ///
		Verbose(integer 0) noLog EForm IRr ///
		EXPosure(varname) OFFset(varname) ///
		D(name) D2 SEParation(string) ///
		GUESS(string) KEEPSINgletons ///
		STANDardize ///
		TAGSEP(name) ZVARname(name) ///
		*]

	// absorb / noabsorb
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

	// Reject absorb() sub-options that ivppmlhdfe cannot honour. The IRLS-IV
	// working variable is not y, so saving FEs / residuals from the inner
	// HDFE solve would produce values for the wrong quantity.
	if strpos(`"`absorb'"', ",") > 0 {
		local _abs_rest = substr(`"`absorb'"', strpos(`"`absorb'"', ",") + 1, .)
		foreach _bad in savefe save residuals resid {
			if strpos(`"`_abs_rest'"', "`_bad'") > 0 {
				di as err "absorb() sub-option `_bad' is not supported by ivppmlhdfe"
				exit 198
			}
		}
	}

	// Detect slope-only absorb (e.g., `absorb(id#c.t)` with no intercept
	// FE). Slope-only absorb is genuinely ambiguous between (a) "include an
	// implicit global constant" and (b) "no constant, slopes only".
	// ppmlhdfe takes (b) and gives a numerically different answer than (a).
	// We refuse cleanly and direct the user to `##c.t` (which both backends
	// agree on bit-for-bit) rather than silently picking one interpretation
	// that disagrees with ppmlhdfe.
	local _abs_pre = `"`absorb'"'
	if strpos(`"`_abs_pre'"', ",") > 0 {
		local _abs_pre = substr(`"`_abs_pre'"', 1, strpos(`"`_abs_pre'"', ",") - 1)
	}
	local _has_intercept_fe = 0
	local _bad_slope_only ""
	foreach _abs_tok of local _abs_pre {
		if regexm("`_abs_tok'", "#c\.") & !regexm("`_abs_tok'", "##c\.") {
			local _bad_slope_only "`_bad_slope_only' `_abs_tok'"
		}
		else {
			local _has_intercept_fe = 1
		}
	}
	if !`_has_intercept_fe' & "`_bad_slope_only'" != "" {
		di as err "slope-only absorb (`_bad_slope_only') not supported by ivppmlhdfe."
		di as err "  ivppmlhdfe and ppmlhdfe disagree on whether slope-only absorb implicitly"
		di as err "  includes a global constant; use {bf:##c.} instead, which yields the same"
		di as err "  fitted coefficient bit-for-bit on both backends."
		exit 198
	}

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

	// d() / d2 option
	// Always materialise the FE-sum into an auto-named tempvar so `predict mu`
	// works without requiring the user to remember `d2` at estimation time.
	// If the user explicitly named the variable, drop any existing copy of it
	// (matching `cap drop` semantics in ppmlhdfe `d2`) so re-fits in loops
	// don't error out with rc=110.
	if "`d2'" != "" & "`d'" == "" {
		local d "_ivppmlhdfe_d"
	}
	if "`d'" == "" {
		local d "_ivppmlhdfe_d"
	}
	capture drop `d'

	// itolerance
	if `itolerance' < 0 {
		local itolerance = -1  // signal to use default in Mata
	}

	// guess()
	if !inlist("`guess'", "", "simple", "mean", "default") {
		di as err "guess() must be simple, mean, or default"
		exit 198
	}

	// tagsep() / zvar()
	if "`tagsep'" != "" {
		confirm new variable `tagsep'
	}
	if "`zvarname'" != "" {
		confirm new variable `zvarname'
	}

	// [*] passthrough: extract display options, keep rest for reghdfe
	_get_diopts diopts options, `options'
	// Remaining options (dof(), pool(), accel(), etc.) will be forwarded
	// to reghdfe via the absorb string using ms_add_comma
	local reghdfe_options `"`options'"'

	// Validate the passthrough against a whitelist of reghdfe options so
	// typos like `tolerence(1e-3)` are rejected with rc=198 instead of being
	// silently forwarded and ignored. ppmlhdfe relies on reghdfe5_parse to
	// reject unknowns, but our Mata code path does not always surface that
	// error cleanly, so we enforce the whitelist here.
	if `"`reghdfe_options'"' != "" {
		local _rh_whitelist "dof pool accel acceleration transform prune keepsingletons precondition"
		// Parse each option: tokens separated by spaces, each is either
		// `name` or `name(args)`. Keep only the name for validation.
		local _rh_tmp `"`reghdfe_options'"'
		local _bad ""
		while `"`_rh_tmp'"' != "" {
			gettoken _opt _rh_tmp : _rh_tmp, parse(" ")
			// Strip trailing `(args)` if present
			local _oname = regexr(`"`_opt'"', "\(.*\)$", "")
			local _oname = strtrim(`"`_oname'"')
			if "`_oname'" == "" continue
			// Match minimum-abbreviation allowed names (up to first match)
			local _ok 0
			foreach _w of local _rh_whitelist {
				if strpos("`_w'", "`_oname'") == 1 {
					local _ok 1
					continue, break
				}
			}
			if !`_ok' local _bad `_bad' `_oname'
		}
		if "`_bad'" != "" {
			di as err "option(s) `_bad' not allowed"
			di as err "  known reghdfe passthrough options: `_rh_whitelist'"
			exit 198
		}
	}

	// separation (aliases match ppmlhdfe)
	local separation : subinstr local separation "ir" "relu", word
	if inlist("`separation'", "", "def", "default", "standard", "on", "auto") {
		local separation "fe simplex relu"
	}
	if inlist("`separation'", "all", "full") {
		local separation "fe simplex relu mu"
	}
	if inlist("`separation'", "no", "off", "none") {
		local separation ""
	}
	// Validate technique names
	local valid_sep "fe simplex relu mu"
	if "`separation'" != "" {
		local bad_sep : list separation - valid_sep
		if "`bad_sep'" != "" {
			di as err "separation(): unknown technique(s) `bad_sep'"
			di as err "  valid: fe simplex relu mu (or aliases: default standard all none)"
			exit 198
		}
	}
	local do_separation = ("`separation'" != "")
	local do_sep_fe      = (strpos("`separation'", "fe") > 0)
	local do_sep_simplex = (strpos("`separation'", "simplex") > 0)
	local do_sep_relu    = (strpos("`separation'", "relu") > 0)
	local do_sep_mu      = (strpos("`separation'", "mu") > 0)

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
		// Accept Stata-standard abbreviations: r/ro/rob/.../robust,
		// cl/clu/.../cluster, un/una/.../unadjusted (mapped to robust since
		// PPML always uses sandwich VCE — matches ppmlhdfe convention).
		// `conventional` and `ols` are also mapped to robust for ppmlhdfe parity.
		if regexm("`vtype'", "^cl") | "`vtype'" == "cluster" local vtype "cluster"
		else if regexm("`vtype'", "^rob") | "`vtype'" == "r" | "`vtype'" == "ro" local vtype "robust"
		else if regexm("`vtype'", "^un") | "`vtype'" == "conventional" | "`vtype'" == "ols" local vtype "robust"
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
	// Support factor variables (i.x) and time-series operators (L.x)
	// via ms_expand_varlist (from ftools, same approach as ppmlhdfe)
	unab depvar : `depvar'

	// Expand exog, endog, instruments through fvexpand + ms_expand_varlist
	// This turns i.industry into 2bn.industry 3bn.industry etc.
	tempvar _touse_fv
	mark `_touse_fv' `if' `in'
	markout `_touse_fv' `depvar'

	if "`exog'" != "" {
		ms_expand_varlist `exog' if `_touse_fv'
		local exog `r(varlist)'
	}
	ms_expand_varlist `endog' if `_touse_fv'
	local endog `r(varlist)'
	ms_expand_varlist `instruments' if `_touse_fv'
	local instruments `r(varlist)'
	if "`clustvars'" != "" {
		local clustvars_unab ""
		local clustvars_display ""
		foreach cv of local clustvars {
			// Handle interaction notation (e.g., exp#imp -> create group var)
			if strpos("`cv'", "#") > 0 {
				tempvar _clust_interact
				local cv_clean : subinstr local cv "#" " ", all
				quietly egen `_clust_interact' = group(`cv_clean')
				local clustvars_unab `clustvars_unab' `_clust_interact'
				local clustvars_display `clustvars_display' `cv'
			}
			else {
				unab cv : `cv'
				// Auto-encode string clusters to a numeric tempvar (matches ppmlhdfe).
				// st_data() in Mata cannot read string variables, so passing a string
				// cluster directly produces "no observations" downstream.
				capture confirm string variable `cv'
				if !_rc {
					tempvar _clust_str
					quietly egen long `_clust_str' = group(`cv')
					local clustvars_unab `clustvars_unab' `_clust_str'
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
		di as err "equation not identified"
		exit 481
	}

	// Reject the empty model: ivppmlhdfe needs at least one regressor (or
	// one absorbed FE that isn't the noabsorb stub). Mata's IRLS-IV builds
	// empty matrices for K=0 and crashes with subscript errors.
	if `n_exog' + `n_endog' == 0 {
		di as err "must specify at least one regressor (an empty model is not supported)"
		exit 198
	}

	// Reject IV that equals an endogenous regressor (or one of the exog).
	// Without this guard the 2SLS first stage degenerates to identity and
	// the estimator silently returns OLS-PPML coefficients.
	if `n_inst' > 0 & `n_endog' > 0 {
		local _bad_iv : list instruments & endog
		if "`_bad_iv'" != "" {
			di as err "instrument cannot equal endogenous regressor: `_bad_iv'"
			exit 198
		}
		local _bad_iv : list instruments & exog
		if "`_bad_iv'" != "" {
			di as err "instrument cannot equal exogenous regressor: `_bad_iv'"
			exit 198
		}
	}

	// Reject constant (zero-variance) instruments — first stage is degenerate.
	if `n_inst' > 0 {
		foreach _iv of local instruments {
			capture confirm numeric variable `_iv'
			if !_rc {
				quietly summarize `_iv', meanonly
				if r(min) == r(max) {
					di as err "instrument `_iv' has zero variance (constant column)"
					exit 459
				}
			}
		}
	}

	// mark sample
	tempvar touse
	mark `touse' `if' `in'
	markout `touse' `depvar' `regressors' `instruments' `clustvars'
	if "`weight'" != "" markout `touse' `exp'

	// Validate exposure>0 BEFORE markout drops missing-offset obs (matches ppmlhdfe).
	// `offset' is the tempvar `offset_var' = ln(`exposure'), which is missing where
	// `exposure' <= 0; if we markout first, the count check would always see zero.
	if "`exposure'" != "" {
		quietly count if `exposure' <= 0 & `touse'
		if r(N) > 0 {
			di as err "exposure() must be greater than zero"
			exit 459
		}
	}
	if "`offset'" != "" markout `touse' `offset'

	quietly count if `depvar' < 0 & `touse'
	if r(N) > 0 {
		di as err "`depvar' must be nonnegative"
		exit 459
	}

	// Reject negative weights — silently propagating them into IRLS leads to
	// contaminated coefficients with no warning. ppmlhdfe rejects with rc=402.
	if "`weight'" != "" {
		quietly count if `exp' < 0 & `touse'
		if r(N) > 0 {
			di as err "negative weights encountered"
			exit 402
		}
		// Drop zero-weight rows from the estimation sample. markout above only
		// drops missing weights, not zeros — but zero-weight rows contribute
		// nothing to score/Hessian, inflate e(N)/e(df) if retained, and can
		// crash Mata paths that expect strictly positive weights. Matches
		// ppmlhdfe's `marksample touse, strok` behavior.
		quietly replace `touse' = 0 if `exp' == 0 & `touse'
	}

	quietly count if `touse'
	local N_full = r(N)
	if `N_full' == 0 error 2000

	// G=1 demotion: if user asked for cluster VCE with a single cluster var
	// that has only 1 group within the touse sample, the cluster sandwich is
	// undefined (df = G - 1 = 0). Demote to robust BEFORE Mata so display,
	// e(df), e(N_clust), and the iter loop all behave consistently.
	if "`vcetype'" == "cluster" & `n_clustvars' == 1 {
		tempvar _g1_check
		quietly egen long `_g1_check' = group(`clustvars') if `touse'
		quietly summarize `_g1_check' if `touse', meanonly
		if r(max) <= 1 {
			di as err "warning: only 1 cluster in `clustvars_display'; switching to vce(robust)"
			local vcetype "robust"
			local clustvars ""
			local clustvars_unab ""
			local clustvars_display ""
			local n_clustvars = 0
		}
	}

	// ---------- 4.  Separation detection ----------
	// FE separation + singletons are handled by reghdfe in Mata via
	// the iweight trick (passing depvar as weight to fixed_effects()).
	// This correctly handles interaction FEs like i.exp#i.imp.
	// Advanced separation (simplex/relu) is also handled in Mata.
	local n_separated = 0
	local n_sep_fe = 0
	local n_sep_advanced = 0

	// ---------- 5.  Singleton detection ----------
	// reghdfe handles both singletons AND FE separation when we pass
	// iweight + depvar to fixed_effects() (see Mata section)

	quietly count if `touse'
	local N_after_sep = r(N)
	if `N_after_sep' == 0 {
		di as err "no observations remaining after separation detection"
		exit 2000
	}

	// ---------- 6.  d() variable ----------
	if "`d'" != "" {
		quietly gen double `d' = .
	}

	// nolog suppresses IRLS iteration output (same as ppmlhdfe)
	if "`log'" == "nolog" & `verbose' >= 0 {
		local verbose = -1
	}

	// ---------- 7.  Display pre-estimation ----------
	if `verbose' > 0 {
		di as txt _n "ivppmlhdfe: IV-PPML with HDFE"
		di as txt "  Dep var:      `depvar'"
		di as txt "  Exogenous:    `exog'"
		di as txt "  Endogenous:   `endog'"
		di as txt "  Instruments:  `instruments'"
		di as txt "  Absorb:       `absorb'"
		di as txt "  VCE:          `vcetype' `clustvars'"
		if "`offset_display'" != "" di as txt "  Offset:       `offset_display'"
		if `n_separated' > 0 {
			di as txt "  Separated:    `n_separated' obs dropped" ///
				" (fe=`n_sep_fe', adv=`n_sep_advanced')"
		}
	}

	// ---------- 8.  tagsep: separation-only mode ----------
	// If tagsep() is specified, run separation detection only (no IRLS).
	// Create HDFE, run FE + simplex + relu, tag separated obs, exit.
	if "`tagsep'" != "" {
		// Run Mata just for HDFE creation + separation
		// We use a lightweight approach: create HDFE, mark separated obs
		tempvar touse_orig
		quietly gen byte `touse_orig' = `touse'

		// Create HDFE with separation (same as full estimation)
		mata: ivppmlhdfe_tagsep( ///
			"`depvar'", "`regressors'", ///
			"`absorb'", "`touse'", "`weight'", "`exp'", ///
			`tolerance', `verbose', ///
			`do_sep_fe', `do_sep_simplex', `do_sep_relu', ///
			"`keepsingletons'" != "", "`zvarname'" ///
		)

		// Create tagsep variable: 1 if in original sample but dropped
		quietly gen byte `tagsep' = `touse_orig' & !`touse'
		label var `tagsep' "Separated observation (ivppmlhdfe)"

		// Report
		quietly count if `tagsep' == 1
		local n_tagged = r(N)
		quietly count if `touse_orig' & !`touse' & `tagsep' == 0
		di as txt _n "Separation detection (no estimation):"
		di as txt "  Tagged as separated: `n_tagged' observations"
		di as txt "  Singletons dropped:  " ivppmlhdfe_tagsep_nsing
		if "`zvarname'" != "" {
			di as txt "  Certificate saved:   `zvarname'"
		}
		exit
	}

	// ---------- 9.  IRLS-IV via Mata ----------
	// Forward reghdfe passthrough options (dof, pool, accel, etc.) via absorb string
	if `"`reghdfe_options'"' != "" {
		ms_add_comma, cmd(`absorb') opt(`reghdfe_options') loc(absorb)
	}

	tempname b_out V_out

	mata: ivppmlhdfe_irls( ///
		"`depvar'", "`exog'", "`endog'", "`instruments'", ///
		"`absorb'", "`touse'", "`weight'", "`exp'", ///
		"`vcetype'", "`clustvars'", ///
		`tolerance', `maxiterations', `verbose', ///
		"`b_out'", "`V_out'", ///
		"`offset'", "`d'", ///
		`do_sep_fe', `do_sep_simplex', `do_sep_relu', ///
		`itolerance', "`guess'", "`keepsingletons'" != "", ///
		`do_sep_mu', "`standardize'" != "" ///
	)

	// Refresh variable lists in case the two-stage collinearity removal
	// inside Mata dropped factor-var redundancies (it rewrites the locals
	// `exog', `endog', `instruments' via st_local). Downstream code uses
	// these for colnames, df_m, Wald chi2, etc. — they must reflect the
	// post-drop varlists.
	local n_exog  : word count `exog'
	local n_endog : word count `endog'
	local n_inst  : word count `instruments'
	local regressors `exog' `endog'

	// ---------- 9.  Post results ----------
	// Loud-failure guard: if IRLS hit maxiter without converging, refuse to
	// post e(b)/e(V) and exit rc=430. Previously we silently posted the
	// non-converged estimate, leaving users to check e(converged) manually
	// — a silent-wrong-answer path. rc=430 matches Stata's standard
	// convention for non-convergent optimizers.
	if ivppmlhdfe_converged == 0 {
		di as err "ivppmlhdfe: IRLS failed to converge in " ivppmlhdfe_iterations " iterations"
		di as err "  coefficients are numerically meaningless and will not be posted."
		di as err "  Try: increasing maxiterations(), loosening tolerance(), simpler FE,"
		di as err "       standardize, or dropping near-separated obs."
		capture mata: mata drop ivppmlhdfe_HDFE
		exit 430
	}

	local N_final = ivppmlhdfe_N
	local n_sep_fe = ivppmlhdfe_num_sep_fe
	local n_sep_advanced = ivppmlhdfe_num_sep_advanced
	local n_sep_mu = ivppmlhdfe_num_sep_mu
	local n_separated = `n_sep_fe' + `n_sep_advanced' + `n_sep_mu'
	local colnames `regressors' _cons
	matrix colnames `b_out' = `colnames'
	matrix colnames `V_out' = `colnames'
	matrix rownames `V_out' = `colnames'

	// Step 1: ereturn post (clears all e())
	
	// Mark omitted/base factor variables (Fix B)
	capture _ms_findomitted `b_out' `V_out'

	ereturn post `b_out' `V_out', esample(`touse') depname(`depvar') buildfvinfo

	// Step 2: Key scalars (before HDFE.post() which may overwrite some)
	ereturn scalar N           = `N_final'
	ereturn scalar converged   = ivppmlhdfe_converged
	ereturn scalar iterations  = ivppmlhdfe_iterations
	ereturn scalar ic          = ivppmlhdfe_iterations

	local K = `n_exog' + `n_endog'
	ereturn scalar df_m        = ivppmlhdfe_actual_rank
	ereturn scalar rank        = ivppmlhdfe_actual_rank

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

	// Cluster df (needed before HDFE.post)
	if "`vcetype'" == "cluster" {
		ereturn scalar df = ivppmlhdfe_G1 - 1
		if `n_clustvars' > 1 {
			local min_G = ivppmlhdfe_G1
			forvalues j = 2/`n_clustvars' {
				local this_G = ivppmlhdfe_G`j'
				if `this_G' < `min_G' local min_G = `this_G'
			}
			ereturn scalar df = `min_G' - 1
		}
		ereturn scalar N_clust = e(df) + 1
	}
	else {
		ereturn scalar df = `N_final' - `K'
	}

	ereturn local cmd          "ivppmlhdfe"

	// Step 3: Call HDFE.post() for metadata (following ppmlhdfe approach)
	// This sets e(vce), e(vcetype), e(clustvar), e(N_clust), e(title2),
	// e(title3), e(depvar), e(indepvars), e(footnote), e(absvars), etc.
	mata: ivppmlhdfe_HDFE.post()

	// HDFE.post() sets e(estat_cmd) = "reghdfe5_estat", which whitelists only
	// reghdfe / ppmlhdfe and refuses to run after ivppmlhdfe (rc=301 on every
	// estat subcommand). Clear it so estat falls back to its default handler.
	ereturn local estat_cmd ""

	// Step 4: Overwrite fields that HDFE.post() set incorrectly for IV-PPML
	// (same pattern as ppmlhdfe.ado lines 330-380)
	ereturn scalar N           = `N_final'
	ereturn scalar ll          = ivppmlhdfe_ll
	ereturn scalar ll_0        = ivppmlhdfe_ll0
	ereturn scalar deviance    = ivppmlhdfe_deviance
	ereturn scalar r2_p        = 1 - ivppmlhdfe_ll / ivppmlhdfe_ll0
	ereturn scalar N_full      = `N_full'
	ereturn scalar num_separated = `n_separated'
	ereturn scalar num_sep_fe = `n_sep_fe'
	ereturn scalar num_sep_advanced = `n_sep_advanced'
	ereturn scalar num_sep_mu = `n_sep_mu'
	ereturn scalar num_singletons = ivppmlhdfe_num_singletons
	ereturn scalar converged   = ivppmlhdfe_converged
	ereturn scalar ic          = ivppmlhdfe_iterations
	ereturn scalar chi2        = e(chi2)

	// Fix J: correct e(df) for robust case using df_a (now available from HDFE.post)
	// Note: don't subtract 1 for _cons — it is absorbed by FE (matching ppmlhdfe)
	if "`vcetype'" != "cluster" {
		ereturn scalar df = `N_final' - `K' - e(df_a)
	}
	else {
		// Re-assert min_G after HDFE.post (matches reghdfe/ppmlhdfe convention)
		local min_G = ivppmlhdfe_G1
		forvalues j = 2/`n_clustvars' {
			local this_G = ivppmlhdfe_G`j'
			if `this_G' < `min_G' local min_G = `this_G'
		}
		ereturn scalar N_clust = `min_G'
		ereturn scalar df = `min_G' - 1
	}

	ereturn local cmd          "ivppmlhdfe"
	ereturn local cmdline      `"ivppmlhdfe `full_cmdline'"'
	ereturn local title        "IV-PPML with High-Dimensional Fixed Effects"
	ereturn local chi2type     "Wald"
	ereturn local predict      "ivppmlhdfe_p"
	ereturn local marginsok    "default"
	ereturn local properties   "b V"
	ereturn local separation   "`separation'"
	ereturn local exogvars     "`exog'"
	ereturn local endogvars    "`endog'"
	ereturn local instruments  "`instruments'"
	if "`offset_display'" != "" ereturn local offset "`offset_display'"
	if "`d'" != "" ereturn local d "`d'"
	// Override cluster variable display names (tempvar -> original names)
	if "`vcetype'" == "cluster" & "`clustvars_display'" != "" {
		ereturn local clustvar "`clustvars_display'"
		local j = 1
		foreach cv of local clustvars_display {
			ereturn local clustvar`j' "`cv'"
			local ++j
		}
	}
	// Fix noabsorb display: override tempvar in e(absvars) with _cons
	if "`noabsorb'" != "" {
		ereturn local absvars "_cons"
		ereturn local title "IV-PPML regression"
		ereturn local title2 ""
	}

	// Blank out OLS-specific fields that HDFE.post() added
	// Keep e(rss), e(rmse), e(ic2) — set by HDFE.post(), matching ppmlhdfe
	ereturn local tss = ""
	ereturn local tss_within = ""
	ereturn local mss = ""
	ereturn local r2 = ""
	ereturn local r2_a = ""
	ereturn local r2_a_within = ""
	ereturn local r2_within = ""
	ereturn local report_constant = ""
	ereturn local sumweights = ""
	ereturn local F = ""
	ereturn local rss = ""

	// Blank out internal IRLS weight leak from HDFE.post()
	// (HDFE stores "aweight" + placeholder for IRLS weights internally;
	// don't expose this as user-facing weight specification)
	if "`weight'" == "" {
		ereturn local wtype = ""
		ereturn local wexp = ""
	}
	else {
		ereturn local wtype "`weight'"
		ereturn local wexp "= `exp'"
	}

	// margins compatibility: set marginsnotok (matching ppmlhdfe)
	ereturn local marginsnotok "stdp Anscombe Cooksd Deviance Hat Likelihood Pearson Response Score Working ADJusted STAndardized STUdentized MODified"

	// Populate e(rss) and e(rmse) for esttab / estat ic compatibility.
	// Must run AFTER `ereturn local predict "ivppmlhdfe_p"` so `predict mu`
	// dispatches to our predict handler (mu is not a default predict option).
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

	// Clean up external HDFE object
	mata: mata drop ivppmlhdfe_HDFE

	// ---------- 10.  Display via Replay ----------
	// nolog suppresses IRLS iteration log (via verbose), NOT the results table
	Replay, `eform' `diopts'
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
// Include ppmlhdfe separation detection (simplex/relu methods)
// ==========================================================================
cap findfile "ppmlhdfe_functions.mata"
if (!_rc) {
    include "`r(fn)'"
    cap findfile "ppmlhdfe_separation_simplex.mata"
    if (!_rc) include "`r(fn)'"
    cap findfile "ppmlhdfe_separation_relu.mata"
    if (!_rc) include "`r(fn)'"
}


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
// Cluster-robust VCE  (returns meat only, without d.f. correction or bread)
// --------------------------------------------------------------------------

real matrix ivppmlhdfe_clustmeat(
	real matrix    scores,
	real colvector clust_id,
	real scalar    K)
{
	real matrix meat, info
	real scalar G, g
	real colvector sg

	info = panelsetup(clust_id, 1)
	G    = rows(info)
	meat = J(K, K, 0)
	for (g = 1; g <= G; g++) {
		sg   = colsum(panelsubmatrix(scores, g, info))'
		meat = meat + sg * sg'
	}
	return(meat)
}


// --------------------------------------------------------------------------
// Create interaction group IDs from two group vectors
// --------------------------------------------------------------------------

real colvector ivppmlhdfe_interact_id(real colvector id1, real colvector id2)
{
	real scalar N
	real colvector srt, result, s1, s2, change

	N = rows(id1)
	srt = order((id1, id2), (1, 2))
	s1 = id1[srt]
	s2 = id2[srt]
	change = 1 \ ((s1[2::N] :!= s1[1::N-1]) :| (s2[2::N] :!= s2[1::N-1]))
	result = J(N, 1, .)
	result[srt] = runningsum(change)
	return(result)
}


// --------------------------------------------------------------------------
// Separation-only mode (for tagsep)
// --------------------------------------------------------------------------

void ivppmlhdfe_tagsep(
	string scalar depvar_s,
	string scalar regressors_s,
	string scalar absorb_s,
	string scalar touse_s,
	string scalar wtype_s,
	string scalar wvar_s,
	real   scalar tolerance,
	real   scalar verbose,
	real   scalar do_sep_fe,
	real   scalar do_sep_simplex,
	real   scalar do_sep_relu,
	real   scalar keepsingletons,
	string scalar zvar_s)
{
	class FixedEffects scalar HDFE
	real scalar drop_sing, num_singletons, K, N
	real colvector y, w_user
	real matrix X
	string scalar absorb_opts
	string rowvector reg_vars

	drop_sing = keepsingletons ? 0 : .
	// Append `precondition` to the absorb option string. If the user already
	// passed reghdfe options via the [*] passthrough, the absorb string is
	// already comma-separated; appending ", precondition" then yields a
	// double-comma that reghdfe5_parse rejects with "invalid 'precondition'".
	// Use a space before precondition iff the absorb string already contains
	// a comma; otherwise insert the comma ourselves.
	if (do_sep_relu) {
		if (strpos(absorb_s, ",")) {
			absorb_opts = absorb_s + " precondition"
		}
		else {
			absorb_opts = absorb_s + ", precondition"
		}
	}
	else {
		absorb_opts = absorb_s
	}

	// Create HDFE with iweight trick for FE separation
	if (do_sep_fe) {
		HDFE = fixed_effects(absorb_opts, touse_s, "iweight", depvar_s, drop_sing, max((0, verbose - 1)))
	}
	else {
		HDFE = fixed_effects(absorb_opts, touse_s, "", "", drop_sing, max((0, verbose - 1)))
	}
	num_singletons = HDFE.num_singletons
	HDFE.save_touse()

	// Load data for simplex/relu
	y = st_data(HDFE.sample, depvar_s)
	N = rows(y)
	reg_vars = tokens(regressors_s)
	K = cols(reg_vars)
	if (K > 0) {
		X = st_data(HDFE.sample, reg_vars)
	}
	else {
		X = J(N, 0, .)
	}
	w_user = (wtype_s != "") ? st_data(HDFE.sample, wvar_s) : J(N, 1, 1)

	// Initialize HDFE weights
	HDFE.load_weights("aweight", "<placeholder>", y, 1)
	HDFE.tolerance = max((1e-4, tolerance))

	// Run simplex/relu if requested and available
	if ((do_sep_simplex | do_sep_relu) & K > 0) {
		if (findexternal("simplex_fix_separation()") != NULL) {
			real rowvector stdev_x
			real colvector non_sep
			real scalar n_drop, target_tol

			stdev_x = diagonal(cholesky(diag(quadvariance(X, w_user))))'
			target_tol = max((1e-12, 0.1 * tolerance))
			HDFE.indepvars = reg_vars

			// Pass scalar 1 (not the length-N vector) when no user weights, so
			// trim_separated_obs / select_not_collinear handle broadcasting.
			real colvector w_for_sep
			w_for_sep = (wtype_s != "") ? w_user : J(1, 1, 1)

			if (do_sep_simplex) {
				non_sep = .
				n_drop = simplex_fix_separation(HDFE, y, X, K, stdev_x,
					w_for_sep, wtype_s != "" ? wtype_s : "", wtype_s != "" ? wvar_s : "",
					target_tol, 1e-12, 1000, non_sep, verbose)
				if (n_drop > 0 & non_sep != .) {
					if (wtype_s != "") w_user = w_user[non_sep]
					else               w_user = J(rows(y), 1, 1)
					N = rows(y)
					stdev_x = diagonal(cholesky(diag(quadvariance(X, w_user))))'
					w_for_sep = (wtype_s != "") ? w_user : J(1, 1, 1)
				}
			}

			if (do_sep_relu) {
				non_sep = .
				n_drop = relu_fix_separation(HDFE, y, X, K, stdev_x,
					w_for_sep, wtype_s != "" ? wtype_s : "", wtype_s != "" ? wvar_s : "",
					target_tol, 1e-4, 1e-8, 100,
					"", zvar_s, 0, 0,
					non_sep, 0, 0, verbose)
				if (n_drop > 0 & non_sep != .) {
					if (wtype_s != "") w_user = w_user[non_sep]
					else               w_user = J(rows(y), 1, 1)
					N = rows(y)
				}
			}
		}
	}

	// Sync touse with final HDFE.sample
	HDFE.save_touse()

	st_numscalar("ivppmlhdfe_tagsep_nsing", num_singletons)
}


// --------------------------------------------------------------------------
// Robust VCE
// --------------------------------------------------------------------------

real matrix ivppmlhdfe_robustvce(
	real matrix    Xhat,
	real matrix    X_dm,
	real colvector w_bread,
	real colvector w_meat,
	real colvector resid,
	real scalar    K,
	real scalar    N)
{
	// w_bread: bread weighting (irls_w = w_user * mu in all cases)
	// w_meat:  per-row score weighting in the sandwich meat
	//          - pw/aw/unweighted: w_meat = irls_w (so meat scales as w_user^2)
	//          - fw:               w_meat = sqrt(w_user) * mu (meat scales as w_user)
	// This matches the ppmlhdfe convention where fweight rows behave as if
	// physically replicated w_user times.
	real matrix bread, meat, V

	bread = invsym(quadcross(Xhat, w_bread, X_dm))
	meat  = quadcross(Xhat :* (w_meat :* resid), Xhat :* (w_meat :* resid))
	V     = (N / (N - 1)) :* bread * meat * bread
	return(V)
}


// --------------------------------------------------------------------------
// Main IRLS-IV loop
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
	| string scalar offset_s,
	string scalar d_var_s,
	real   scalar do_sep_fe,
	real   scalar do_sep_simplex,
	real   scalar do_sep_relu,
	real   scalar user_itolerance,
	string scalar guess_s,
	real   scalar keepsingletons,
	real   scalar do_sep_mu,
	real   scalar do_standardize)
{
	// ---- Declarations ----
	real colvector y, mu, eta, old_eta, z, irls_w, resid, w_user
	real colvector offset
	real matrix X, Z_excl, data
	real matrix X_dm, Z_dm, Xhat
	real colvector z_dm, b, b_full
	real matrix V, V_slope
	real scalar N, K_exog, K_endo, K, L, K_total, n_inst
	real scalar iter, converged, ok
	real scalar deviance, old_deviance, delta_dev, eps, denom_eps
	real colvector b_old
	real scalar beta_change
	real scalar mean_y, b_cons
	string rowvector exog_vars, endog_vars, inst_vars
	real scalar has_cluster, has_weight, has_offset
	real colvector sort_order
	real matrix ZwZ_f, Pi_f

	// Step-halving
	real scalar iter_step_halving, num_step_halving, max_step_halving
	real scalar step_halving_memory

	// Adaptive tolerance
	real scalar start_inner_tol, target_inner_tol, alt_tol

	// HDFE object (from reghdfe)
	class FixedEffects scalar HDFE

	// Multi-way clustering
	string rowvector clust_varnames
	real scalar n_clust
	real matrix clust_ids

	// Standardization
	real rowvector stdev_x, stdev_z
	real scalar stdev_y

	// Log-likelihood
	real scalar ll, ll_0, ll_0_mu

	// Singletons
	real scalar num_singletons

	// Separation detection
	real scalar num_sep_fe, num_sep_advanced, num_sep_mu
	real scalar log_septol, adjusted_log_septol, min_eta_pos
	real colvector sep_mask

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
	has_offset  = (args() >= 16 & offset_s != "")

	// Separation flags (optional args with defaults)
	if (args() < 19) do_sep_fe = 1
	if (args() < 20) do_sep_simplex = 0
	if (args() < 21) do_sep_relu = 0
	if (args() < 22) user_itolerance = -1
	if (args() < 23) guess_s = ""
	if (args() < 24) keepsingletons = 0
	if (args() < 25) do_sep_mu = 0
	if (args() < 26) do_standardize = 0
	num_sep_fe = 0
	num_sep_advanced = 0
	num_sep_mu = 0

	// keepsingletons warning (matches ppmlhdfe.mata:197 convention).
	if (keepsingletons & verbose > -1) {
		printf("{err}warning: keeping singleton groups will keep fixed effects that cause separation\n")
	}

	// ---- Parse cluster variables ----
	if (has_cluster) {
		clust_varnames = tokens(clustvar_s)
		n_clust = cols(clust_varnames)
	}
	else {
		n_clust = 0
	}

	// Step-halving parameters. Match ppmlhdfe.mata defaults
	// (max_step_halving=2, step_halving_memory=0.9).
	max_step_halving = 2
	step_halving_memory = 0.9

	// Adaptive tolerance parameters
	start_inner_tol = 1e-4
	target_inner_tol = user_itolerance > 0 ? user_itolerance : max((1e-12, 0.1 * tolerance))

	// ---- Create HDFE object ----
	// When FE separation is on, use the iweight trick (following ppmlhdfe):
	// pass depvar as iweight to fixed_effects(). reghdfe then drops FE groups
	// where sum(weight)==0, i.e., where all y==0. This correctly handles
	// interaction FEs like i.exp#i.imp that the old manual egen approach missed.
	{
		real scalar N_before_hdfe, drop_sing
		string scalar absorb_opts
		N_before_hdfe = sum(st_data(., touse_s))
		drop_sing = keepsingletons ? 0 : .
		// Append precondition when relu is requested (needed by solve_lse).
		// If user passed reghdfe options via [*] passthrough, absorb_s already
		// contains a comma; appending ", precondition" would produce a double
		// comma that reghdfe5_parse rejects. See tagsep branch above for the
		// same fix.
		if (do_sep_relu) {
			if (strpos(absorb_s, ",")) {
				absorb_opts = absorb_s + " precondition"
			}
			else {
				absorb_opts = absorb_s + ", precondition"
			}
		}
		else {
			absorb_opts = absorb_s
		}
		if (do_sep_fe) {
			HDFE = fixed_effects(absorb_opts, touse_s, "iweight", depvar_s, drop_sing, max((0, verbose - 1)))
		}
		else {
			HDFE = fixed_effects(absorb_opts, touse_s, "", "", drop_sing, max((0, verbose - 1)))
		}
		num_singletons = HDFE.num_singletons
		// Count FE-separated obs (total dropped minus singletons)
		num_sep_fe = N_before_hdfe - rows(HDFE.sample) - num_singletons
		if (num_sep_fe < 0) num_sep_fe = 0
		// Sync Stata touse with HDFE.sample (needed for d() and e(sample))
		HDFE.save_touse()
	}

	// ---- Load data using HDFE.sample (guaranteed consistent with HDFE internals) ----
	// NOTE: Must use HDFE.sample, not touse_s, because fixed_effects()
	// may drop singletons that are reflected in HDFE.sample but not
	// perfectly synchronized with touse_s. This follows ppmlhdfe's approach.
	y = st_data(HDFE.sample, depvar_s)
	N = rows(y)

	if (K_exog > 0 & K_endo > 0) {
		X = st_data(HDFE.sample, (exog_vars, endog_vars))
	}
	else if (K_exog > 0) {
		X = st_data(HDFE.sample, exog_vars)
	}
	else if (K_endo > 0) {
		X = st_data(HDFE.sample, endog_vars)
	}
	else {
		X = J(N, 0, .)
	}

	if (n_inst > 0) {
		Z_excl = st_data(HDFE.sample, inst_vars)
	}
	else {
		Z_excl = J(N, 0, .)
	}

	if (has_weight) {
		w_user = st_data(HDFE.sample, wvar_s)
	}
	else {
		w_user = J(N, 1, 1)
	}

	// Load cluster IDs
	if (has_cluster) {
		clust_ids = st_data(HDFE.sample, clust_varnames)
	}

	// Load offset. We center it to weighted mean 0 for the IRLS loop — this
	// keeps the initial mu finite when the user passes a large absolute-scale
	// offset (e.g., log(trade_$) ≈ 25), which would otherwise make
	// mu_init = exp(offset) ≈ 1e11 and blow up the bread matrix.  ppmlhdfe
	// avoids this via compute_constant=1 inside partial_out; we keep
	// compute_constant=0, so we re-center manually.  After IRLS converges,
	// `eta` and `mu` already equal log(mu_user) and mu_user respectively
	// (they absorbed the centring shift into the fitted b_cons_c) — so we
	// only need to restore `offset` to its original values for the
	// post-processing formulas (b_cons, d_vals, ll_0).
	real scalar offset_mean
	real colvector offset_orig
	offset_mean = 0
	if (has_offset) {
		offset_orig = st_data(HDFE.sample, offset_s)
		offset_mean = quadsum(w_user :* offset_orig) / quadsum(w_user)
		offset = offset_orig :- offset_mean
	}
	else {
		offset_orig = J(N, 1, 0)
		offset = J(N, 1, 0)
	}

	// ---- Initialize HDFE weights ----
	// Use y as placeholder weight vector (following ppmlhdfe pattern).
	HDFE.load_weights("aweight", "<placeholder for mu>", y, 1)
	HDFE.tolerance = max((start_inner_tol, tolerance))

	// ---- Two-stage collinearity removal (mirrors ppmlhdfe pattern) ----
	// Detect linearly dependent columns within exog, endog, and instrument
	// blocks AFTER projecting them through the FE structure. ppmlhdfe's
	// remove_collinears (ppmlhdfe_functions.mata:251) runs once on the exog
	// block; we run it three times — once per IV block — because IV-PPML has
	// three independent matrices that must each be full rank for the 2SLS
	// solve to be well-posed (X'WX bread, Z'WZ first-stage cross product).
	//
	// Without this step, factor-variable expansions like i.x + absorb(x),
	// i.x##i.y, ibn.x, etc. silently produce zero/redundant columns that
	// blow up the IRLS bread invsym (line 1845) or trip the cascade
	// guardrail (rc=9003). Surviving column names are written back to
	// Stata locals so e(b) colnames stay aligned.
	// NOTE: no findexternal() gate — that returns NULL for functions
	// included inside .ado scope (they live in the .ado's local Mata
	// namespace, not the global external one). select_not_collinear is
	// reliably available because ppmlhdfe_functions.mata is included at
	// .ado load time (lines 922-933 above). If that include ever fails,
	// the call below will throw a Mata runtime error rather than silently
	// skipping collinearity removal.
	if ((K + n_inst) > 0) {
		real matrix data_for_coll
		real rowvector ok_block, exog_keep, endog_keep, inst_keep
		real colvector w_for_coll
		real scalar j_coll
		real scalar n_dropped_exog, n_dropped_endog, n_dropped_inst

		// Build a temporary copy [X, Z_excl] for partialling out
		if (K > 0 & n_inst > 0) {
			data_for_coll = (X, Z_excl)
		}
		else if (K > 0) {
			data_for_coll = X
		}
		else {
			data_for_coll = Z_excl
		}

		// Partial out the temp copy against the FE structure with USER
		// weights (NOT the IRLS placeholder weights, which haven't been
		// updated yet). Backup HDFE state, swap weights, partial out,
		// restore — same dance as ppmlhdfe_functions.mata:275-289.
		real scalar bk_tol_coll
		string scalar bk_wt_type_coll, bk_wt_var_coll
		real colvector bk_w_coll
		bk_tol_coll = HDFE.tolerance
		bk_wt_type_coll = HDFE.weight_type
		bk_wt_var_coll = HDFE.weight_var
		bk_w_coll = HDFE.weight
		// Always use aweight placeholder + a length-N vector — this matches
		// the existing IRLS pattern at line 1392 and avoids leaving HDFE in
		// a half-initialized state when no user weights are present
		// (HDFE.load_weights("", ...) leaves HDFE.weight as a 1×1 scalar
		// which corrupts subsequent partial_out invocations).
		HDFE.load_weights("aweight", "<placeholder for collinearity>",
		                  has_weight ? w_user : J(N, 1, 1), 0)
		HDFE.tolerance = target_inner_tol
		(void) --HDFE.verbose
		HDFE._partial_out(data_for_coll, 0, 0, 0, 1)
		(void) ++HDFE.verbose
		HDFE.load_weights(bk_wt_type_coll, bk_wt_var_coll, bk_w_coll, 0)
		HDFE.tolerance = bk_tol_coll
		bk_w_coll = .

		w_for_coll = has_weight ? w_user : J(1, 1, 1)

		// ---- Block 1: exog (within-block collinearity) ----
		exog_keep = J(1, K_exog, 1)
		if (K_exog > 0) {
			ok_block = .
			(void) select_not_collinear(data_for_coll[., 1..K_exog],
			                            w_for_coll, ok_block)
			if (ok_block != .) {
				exog_keep = J(1, K_exog, 0)
				for (j_coll = 1; j_coll <= cols(ok_block); j_coll++) {
					exog_keep[ok_block[j_coll]] = 1
				}
			}
		}

		// ---- Block 2: endog (within-block collinearity) ----
		endog_keep = J(1, K_endo, 1)
		if (K_endo > 0) {
			ok_block = .
			(void) select_not_collinear(
				data_for_coll[., (K_exog + 1)..K], w_for_coll, ok_block)
			if (ok_block != .) {
				endog_keep = J(1, K_endo, 0)
				for (j_coll = 1; j_coll <= cols(ok_block); j_coll++) {
					endog_keep[ok_block[j_coll]] = 1
				}
			}
		}

		// ---- Block 3: instruments (Z_excl alone, within-block) ----
		inst_keep = J(1, n_inst, 1)
		if (n_inst > 0) {
			ok_block = .
			(void) select_not_collinear(
				data_for_coll[., (K + 1)..(K + n_inst)],
				w_for_coll, ok_block)
			if (ok_block != .) {
				inst_keep = J(1, n_inst, 0)
				for (j_coll = 1; j_coll <= cols(ok_block); j_coll++) {
					inst_keep[ok_block[j_coll]] = 1
				}
			}
		}

		// ---- Apply drops to original X and Z_excl matrices ----
		n_dropped_exog  = K_exog - sum(exog_keep)
		n_dropped_endog = K_endo - sum(endog_keep)
		n_dropped_inst  = n_inst - sum(inst_keep)

		if (n_dropped_exog + n_dropped_endog + n_dropped_inst > 0) {
			real rowvector exog_idx_c, endog_idx_c, inst_idx_c
			real matrix X_kept, Z_kept
			string rowvector new_exog_names, new_endog_names, new_inst_names

			if (K_exog > 0) {
				exog_idx_c = selectindex(exog_keep)
				new_exog_names = exog_vars[exog_idx_c]
			}
			else {
				new_exog_names = J(1, 0, "")
			}
			if (K_endo > 0) {
				endog_idx_c = selectindex(endog_keep)
				new_endog_names = endog_vars[endog_idx_c]
			}
			else {
				new_endog_names = J(1, 0, "")
			}
			if (n_inst > 0) {
				inst_idx_c = selectindex(inst_keep)
				new_inst_names = inst_vars[inst_idx_c]
			}
			else {
				new_inst_names = J(1, 0, "")
			}

			// Subset the raw (un-partialled) X and Z_excl
			if (cols(new_exog_names) > 0 & cols(new_endog_names) > 0) {
				X_kept = (X[., exog_idx_c],
				          X[., (K_exog :+ endog_idx_c)])
			}
			else if (cols(new_exog_names) > 0) {
				X_kept = X[., exog_idx_c]
			}
			else if (cols(new_endog_names) > 0) {
				X_kept = X[., (K_exog :+ endog_idx_c)]
			}
			else {
				X_kept = J(N, 0, .)
			}
			if (cols(new_inst_names) > 0) {
				Z_kept = Z_excl[., inst_idx_c]
			}
			else {
				Z_kept = J(N, 0, .)
			}

			X = X_kept
			Z_excl = Z_kept
			K_exog = cols(new_exog_names)
			K_endo = cols(new_endog_names)
			K = K_exog + K_endo
			n_inst = cols(new_inst_names)
			L = K_exog + n_inst

			// Update Mata token vectors so downstream separation/IRLS
			// references to exog_vars/endog_vars/inst_vars stay valid
			exog_vars  = new_exog_names
			endog_vars = new_endog_names
			inst_vars  = new_inst_names

			// Identification check after dropping
			if (n_inst < K_endo) {
				printf("{err}error: equation not identified after collinearity removal ")
				printf("{err}(n_inst=%g < K_endo=%g)\n", n_inst, K_endo)
				exit(481)
			}

			// Write surviving names back to Stata so e(b) colnames stay aligned
			st_local("exog",        invtokens(new_exog_names))
			st_local("endog",       invtokens(new_endog_names))
			st_local("instruments", invtokens(new_inst_names))

			if (verbose > -1) {
				if (n_dropped_exog > 0) {
					printf("{txt}note: %g exog variable%s omitted because of collinearity\n",
						n_dropped_exog, n_dropped_exog > 1 ? "s" : "")
				}
				if (n_dropped_endog > 0) {
					printf("{txt}note: %g endog variable%s omitted because of collinearity\n",
						n_dropped_endog, n_dropped_endog > 1 ? "s" : "")
				}
				if (n_dropped_inst > 0) {
					printf("{txt}note: %g instrument%s omitted because of collinearity\n",
						n_dropped_inst, n_dropped_inst > 1 ? "s" : "")
				}
			}
		}

		// Free temporary
		data_for_coll = J(0, 0, .)
	}

	// ---- Advanced separation detection (simplex/relu) ----
	// Uses ppmlhdfe's compiled Mata functions (included at load time)
	// Separation depends on (y, X, FE), NOT on instruments Z
	if ((do_sep_simplex | do_sep_relu) & K > 0) {
		real scalar n_drop_sep, has_ppml_sep, N_before_sep
		real colvector non_sep

		// Check if ppmlhdfe functions are actually available
		has_ppml_sep = (findexternal("simplex_fix_separation()") != NULL)

		if (!has_ppml_sep) {
			printf("{txt}(note: simplex/relu separation requires ppmlhdfe; only fe separation will run)\n")
		}
		else {
			// Compute stdev_x (needed by separation functions)
			stdev_x = diagonal(cholesky(diag(quadvariance(X, w_user))))'

			// Set HDFE.indepvars (required by simplex/relu assertions)
			if (K_exog > 0 & K_endo > 0) {
				HDFE.indepvars = (exog_vars, endog_vars)
			}
			else if (K_endo > 0) {
				HDFE.indepvars = endog_vars
			}
			else if (K_exog > 0) {
				HDFE.indepvars = exog_vars
			}

			// IMPORTANT: when there are no user weights, ppmlhdfe passes a SCALAR
			// J(1,1,1) as `true_w` to the separation routines. trim_separated_obs
			// only trims true_w when weight_var != "", and select_not_collinear
			// calls quadcross(x, w, x) which broadcasts a scalar w. Passing the
			// length-N vector J(N,1,1) instead leaves true_w out of sync with the
			// trimmed x and triggers a 3200 conformability error in quadcross
			// (deepest manifestation: keepsingletons + relu).
			real colvector w_for_sep
			w_for_sep = has_weight ? w_user : J(1, 1, 1)

			// Simplex method: detects separation due to regressors.
			// simplex_fix_separation / relu_fix_separation INTERNALLY call
			// trim_separated_obs (ppmlhdfe_functions.mata:179) which mutates
			// HDFE.sample, calls HDFE.reload(0), handles the cascade where
			// reload detects new singletons, and trims y, X in place. BUT
			// `HDFE = HDFE.reload(0)` inside trim_separated_obs is a LOCAL
			// reassignment that doesn't propagate back through Mata's call
			// stack. Our caller's HDFE has the mutated sample but stale
			// internal factor structure. We explicitly reload HERE so
			// downstream HDFE._partial_out works on a coherent object.
			if (do_sep_simplex) {
				N_before_sep = rows(y)
				non_sep = .
				n_drop_sep = simplex_fix_separation(HDFE, y, X, K, stdev_x,
					w_for_sep,
					has_weight ? wtype_s : "",
					has_weight ? wvar_s : "",
					target_inner_tol, 1e-12, 1000, non_sep, verbose)
				if (n_drop_sep > 0 & non_sep != .) {
					Z_excl = Z_excl[non_sep, .]
					if (has_weight) w_user = w_user[non_sep]
					else            w_user = J(rows(y), 1, 1)
					if (has_offset) {
						offset      = offset[non_sep]
						offset_orig = offset_orig[non_sep]
					}
					else {
						offset      = J(rows(y), 1, 0)
						offset_orig = J(rows(y), 1, 0)
					}
					if (has_cluster) clust_ids = clust_ids[non_sep, .]
					N = rows(y)
					num_sep_advanced = num_sep_advanced + (N_before_sep - N)
					// Reload HDFE so its factors / N match our trimmed sample
					real scalar _ns_s
					_ns_s = HDFE.num_singletons
					HDFE = HDFE.reload(0)
					HDFE.num_singletons = _ns_s
					HDFE.save_touse()  // sync Stata touse for e(sample) + predict
					stdev_x = diagonal(cholesky(diag(quadvariance(X, w_user))))'
					w_for_sep = has_weight ? w_user : J(1, 1, 1)
				}
			}

			// ReLU method: detects separation due to regressors + FE
			if (do_sep_relu) {
				N_before_sep = rows(y)
				non_sep = .
				n_drop_sep = relu_fix_separation(HDFE, y, X, K, stdev_x,
					w_for_sep,
					has_weight ? wtype_s : "",
					has_weight ? wvar_s : "",
					target_inner_tol,
					1e-4, 1e-8, 100,
					"", "", 0, 0,
					non_sep, 0, 0, verbose)
				if (n_drop_sep > 0 & non_sep != .) {
					Z_excl = Z_excl[non_sep, .]
					if (has_weight) w_user = w_user[non_sep]
					else            w_user = J(rows(y), 1, 1)
					if (has_offset) {
						offset      = offset[non_sep]
						offset_orig = offset_orig[non_sep]
					}
					else {
						offset      = J(rows(y), 1, 0)
						offset_orig = J(rows(y), 1, 0)
					}
					if (has_cluster) clust_ids = clust_ids[non_sep, .]
					N = rows(y)
					num_sep_advanced = num_sep_advanced + (N_before_sep - N)
					real scalar _ns_r
					_ns_r = HDFE.num_singletons
					HDFE = HDFE.reload(0)
					HDFE.num_singletons = _ns_r
					HDFE.save_touse()  // sync Stata touse
				}
			}

			if (num_sep_advanced > 0) {
				// Reload HDFE weights with updated y
				HDFE.load_weights("aweight", "<placeholder for mu>", y, 1)
				HDFE.tolerance = max((start_inner_tol, tolerance))
			}
		}
	}

	// ---- Standardize data (optional) ----
	// Divides X and Z by column stdev for numerical stability.
	// Coefficients and VCE are back-transformed after convergence.
	if (do_standardize & K > 0) {
		stdev_x = J(1, K, 1)
		real scalar _kk_std
		for (_kk_std = 1; _kk_std <= K; _kk_std++) {
			stdev_x[_kk_std] = sqrt(quadvariance(X[., _kk_std], w_user))
			if (stdev_x[_kk_std] > 0) {
				X[., _kk_std] = X[., _kk_std] :/ stdev_x[_kk_std]
			}
			else {
				stdev_x[_kk_std] = 1  // constant column: leave as-is
			}
		}
		stdev_z = J(1, n_inst, 1)
		for (_kk_std = 1; _kk_std <= n_inst; _kk_std++) {
			stdev_z[_kk_std] = sqrt(quadvariance(Z_excl[., _kk_std], w_user))
			if (stdev_z[_kk_std] > 0) {
				Z_excl[., _kk_std] = Z_excl[., _kk_std] :/ stdev_z[_kk_std]
			}
			else {
				stdev_z[_kk_std] = 1
			}
		}
		if (verbose > 0) {
			printf("{txt}(data standardized: stdev_x = ")
			for (_kk_std = 1; _kk_std <= K; _kk_std++) printf("%9.4f ", stdev_x[_kk_std])
			printf(")\n")
		}
	}
	else {
		stdev_x = J(1, K, 1)
		stdev_z = J(1, n_inst, 1)
	}
	stdev_y = 1  // y standardization not implemented (Poisson scale-specific)

	// ---- Initialise mu ----
	mean_y = mean(y, w_user)
	if (guess_s == "mean") {
		// guess(mean): start all mu at mean(y)
		mu = J(N, 1, mean_y)
	}
	else {
		// Default / guess(simple): mu = 0.5 * (y + mean(y))
		mu = 0.5 :* (y :+ mean_y)
	}
	// Y-dependent censoring (matching ppmlhdfe's censor_mu):
	// floor at max(0.05*y, 1e-3), not flat 1e-4
	mu = rowmax((mu, rowmax((0.05 :* y, J(N, 1, 1e-3)))))
	eta = log(mu) + offset

	// ---- IRLS-IV loop ----
	converged = 0
	ok = 0
	deviance = .
	eps = .
	beta_change = .
	b     = J(K, 1, 0)
	b_old = J(K, 1, 0)
	iter_step_halving = 0
	num_step_halving = 0
	alt_tol = start_inner_tol
	// Persistent mu-separation mask (accumulated across iterations)
	if (do_sep_mu) sep_mask = J(N, 1, 0)

	if (verbose > -1) {
		printf("{txt}\nIRLS-IV iterations (N = %g, K = %g, L = %g)\n", N, K, L)
		printf("{txt}{hline 60}\n")
	}

	for (iter = 1; iter <= maxiter; iter++) {

		// (a0) Runaway-divergence guard (matching ppmlhdfe.mata line 620 +
		// additional |b|/|eta| sanity bounds to catch the Class C pathology
		// where mu stays finite but IRLS thrashes to garbage coefficients).
		if (hasmissing(mu) | missing(quadsum(mu))) {
			printf("{err}error: mu has infinite or missing values on iteration %g; aborting\n", iter)
			printf("{err}       this typically indicates IRLS divergence on a heavily\n")
			printf("{err}       separated or near-collinear model. Try:\n")
			printf("{err}         - separation(all)        (enables in-loop mu-separation)\n")
			printf("{err}         - standardize            (rescales X, Z columns)\n")
			printf("{err}         - dropping the most fragile FE dimension\n")
			exit(9003)
		}
		if (iter > 10 & K > 0) {
			real scalar _max_b, _max_eta
			_max_b   = max(abs(b))
			_max_eta = max(abs(eta))
			// If |b| or |eta| has exploded beyond any plausible economic
			// magnitude, the IRLS has clearly diverged — bail out before
			// posting nonsense. Thresholds chosen generously to avoid false
			// positives: exp(30)=1e13, |b|=1e6 on standardized data is absurd.
			if (_max_b > 1e6 | _max_eta > 30) {
				printf("{err}error: IRLS diverged at iter %g (max|b|=%9.2e, max|eta|=%9.2e); aborting\n", iter, _max_b, _max_eta)
				printf("{err}       coefficients have grown past sanity bounds. Try:\n")
				printf("{err}         - separation(all)        (in-loop mu-separation)\n")
				printf("{err}         - standardize            (rescale X, Z)\n")
				printf("{err}         - dropping the most fragile FE dimension\n")
				printf("{err}         - a longer panel (Class C needs large N_c, T)\n")
				exit(9003)
			}
		}

		// (a) Working depvar: z = eta - offset - 1 + y/mu
		//     (subtract offset before demeaning, following ppmlhdfe)
		z = eta - offset :- 1 + y :/ mu
		// For y=0 obs: set z = eta - offset - 1 exactly (avoid 0/epsilon noise,
		// matching ppmlhdfe's zero_sample handling during separation check)
		if (do_sep_mu) {
			real colvector zero_idx
			zero_idx = selectindex(y :== 0)
			if (rows(zero_idx) > 0) {
				z[zero_idx] = eta[zero_idx] - offset[zero_idx] :- 1
			}
		}

		// (b) IRLS weights: w = w_user * mu  (floor to prevent HDFE assertion crash)
		irls_w = w_user :* mu
		irls_w = rowmax((irls_w, J(N, 1, 1e-20)))

		// (c) Update HDFE weights
		HDFE.update_sorted_weights(irls_w)
		HDFE.update_cvar_objects()

		// (d) Stack data matrix (z without offset, X, Z_excl)
		data = (z, X, Z_excl)

		// (e) Demean all variables via HDFE
		_edittozerotol(data, min((tolerance, 1e-12)))
		HDFE._partial_out(data, 0, 0, 0, 1)
		// Round near-zero demeaned values to exact zero (prevents
		// floating-point noise accumulation with high-dimensional FE)
		_edittozerotol(data, min((tolerance, 1e-12)))

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
		if (K > 0) b_old = b
		ivppmlhdfe_gmm(z_dm, X_dm, Z_dm, irls_w, b, resid)
		if (K > 0) beta_change = max(abs(b - b_old))

		// (h) Update eta = z - resid + offset
		if (!iter_step_halving) swap(old_eta, eta)
		eta = z - resid + offset

		// (h2) Mu-separation check (matching ppmlhdfe)
		// Detect y=0 obs where eta is diverging to -infinity.
		// Accumulates across iterations (once separated, always separated).
		if (do_sep_mu & iter > 1) {
			real scalar n_new_sep
			log_septol = log(1e-6)  // = -13.82
			if (any(y :> 0)) {
				min_eta_pos = min(select(eta, y :> 0))
				adjusted_log_septol = log_septol + min((min_eta_pos + 5, 0))
			}
			else {
				adjusted_log_septol = log_septol
			}
			// Accumulate: |= so once separated, always separated
			sep_mask = sep_mask :| ((eta :<= adjusted_log_septol) :& (y :== 0))
			n_new_sep = sum(sep_mask) - num_sep_mu
			if (n_new_sep > 0) {
				num_sep_mu = sum(sep_mask)
				if (verbose > -1) {
					printf("{txt}(mu-separation: %g obs detected at iter %g)\n", n_new_sep, iter)
				}
			}
		}

		// (i) Update mu (floor at epsilon(100) ~1.42e-14, matching ppmlhdfe)
		mu = exp(eta)
		// Set mu=0 for separated obs (matching ppmlhdfe), then floor
		if (do_sep_mu & num_sep_mu > 0) {
			real colvector sep_idx
			sep_idx = selectindex(sep_mask)
			mu[sep_idx] = J(num_sep_mu, 1, 0)
		}
		_vector_scalar_max(mu, epsilon(100))

		// (j) Deviance (with numerical safeguards matching ppmlhdfe)
		old_deviance = deviance
		deviance = quadsum((mu - y) :* w_user) ///
			+ quadcross(y, (y :> 0) :* w_user, log(y) - eta)
		if (2 * deviance / N < epsilon(1)) deviance = 0
		deviance = 2 * edittozerotol(deviance, epsilon(1))
		if (deviance < 0) deviance = 0

		// (k) Convergence + step-halving + adaptive tolerance
		if (iter > 1) {
			delta_dev = old_deviance - deviance
			// Clip delta_dev: deviance near zero can't decrease by more than itself
			if (deviance < 0.1 * delta_dev) delta_dev = deviance
			denom_eps = max((min((deviance, old_deviance)), 0.1))
			eps = abs(delta_dev) / denom_eps

			if (eps < tolerance) {
				if (HDFE.tolerance <= 1.1 * target_inner_tol | HDFE.G == 1) {
					ok = ok + 1
					if (ok >= 1) converged = 1
				}
			}
			else if (delta_dev < 0 & num_step_halving < max_step_halving) {
				eta = step_halving_memory * old_eta + (1 - step_halving_memory) * eta
				// Clip eta at -10 on repeated step-halving (matching ppmlhdfe)
				if (num_step_halving > 0) {
					eta = rowmax((eta, J(N, 1, -10)))
				}
				mu = exp(eta)
				_vector_scalar_max(mu, epsilon(100))
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
			if (K > 0 & beta_change < .) printf("{txt}  db = {res}%-9.4e", beta_change)
			printf("{txt}  tol = {res}%5.0e", HDFE.tolerance)
			if (iter_step_halving) printf("{txt}  H")
			if (ok) printf("{txt}  O")
			printf("\n")
		}

		if (iter_step_halving) {
			deviance = old_deviance
			continue
		}

		if (converged) break

		// Adaptive inner tolerance
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

	// Restore offset to the user's original scale for post-processing.
	// eta / mu / irls_w are already on the user's scale (the IRLS absorbed
	// the offset-centring shift into the fitted b_cons_c), so they need no
	// adjustment. See the offset-centering comment at the top of the function.
	offset = offset_orig

	b_cons = mean(eta - offset, irls_w) - mean(X, irls_w) * b
	b_full = b \ b_cons
	K_total = K + 1

	// Recompute Xhat for VCE
	ZwZ_f = invsym(cross(Z_dm, irls_w, Z_dm))
	Pi_f  = ZwZ_f * cross(Z_dm, irls_w, X_dm)
	Xhat  = Z_dm * Pi_f

	// Rank check: detect collinearity via diag0cnt on bread
	real scalar actual_rank
	actual_rank = K

	// Effective sample size: with fweight, observations are replicated w_i
	// times, so the VCE dof_adj and reported N must use sum(w_user). For
	// pweight / aweight / unweighted, effective N is the physical row count.
	real scalar N_eff
	N_eff = (wtype_s == "fweight") ? quadsum(w_user) : N

	// Per-row meat weight for the sandwich VCE.  fweight rows behave as if
	// physically replicated w_user times, so the meat scales as w_user (one
	// power), not w_user^2.  For pw/aw/unweighted the meat scales as the
	// usual w_user^2 (sampling-variance convention).  See robustvce comment.
	real colvector w_meat
	if (wtype_s == "fweight") {
		w_meat = sqrt(w_user) :* mu
	}
	else {
		w_meat = irls_w
	}

	// ---- Compute VCE ----
	if (K == 0) {
		// No regressors: empty VCE
		V_slope = J(0, 0, 0)
	}
	else if (vcetype_s == "cluster") {
		// Compute scores = Xhat * w * resid
		{
			real matrix scores, bread
			real scalar j, mask, subset_size, first_j, sign
			real colvector cid
			real matrix meat_sub, info_sub
			real scalar G_sub, g_sub
			real colvector sg_sub, srt_sub
			real scalar N_clust_min

			bread  = invsym(cross(Xhat, irls_w, X_dm))
			actual_rank = K - diag0cnt(bread)
			scores = Xhat :* (irls_w :* resid)

			if (n_clust == 1) {
				// Single cluster
				srt_sub = order(clust_ids, 1)
				meat_sub = ivppmlhdfe_clustmeat(scores[srt_sub, .], clust_ids[srt_sub], K)
				G_sub = rows(panelsetup(clust_ids[srt_sub], 1))
				if (G_sub <= 1) {
					printf("{err}Warning: only %g cluster(s); cluster VCE not defined. Using robust.\n", G_sub)
					V_slope = ivppmlhdfe_robustvce(Xhat, X_dm, irls_w, w_meat, resid, K, N_eff)
				}
				else {
					V_slope = (G_sub / (G_sub - 1)) :* bread * meat_sub * bread
				}
				st_numscalar("ivppmlhdfe_G1", G_sub)
			}
			else {
				// Multi-way CGM clustering
				V_slope = J(K, K, 0)
				N_clust_min = N

				// Iterate over all non-empty subsets of {1..n_clust}
				for (mask = 1; mask < 2^n_clust; mask++) {
					// Build combined cluster ID for this subset
					subset_size = 0
					first_j = 1
					for (j = 1; j <= n_clust; j++) {
						if (mod(trunc(mask / 2^(j-1)), 2) == 1) {
							subset_size++
							if (first_j) {
								cid = clust_ids[., j]
								first_j = 0
							}
							else {
								cid = ivppmlhdfe_interact_id(cid, clust_ids[., j])
							}
						}
					}

					// Compute meat for this cluster
					srt_sub = order(cid, 1)
					meat_sub = ivppmlhdfe_clustmeat(scores[srt_sub, .], cid[srt_sub], K)
					G_sub = rows(panelsetup(cid[srt_sub], 1))

					// CGM sign: + for odd subset size, - for even
					sign = mod(subset_size, 2) == 1 ? 1 : -1

					// Add with d.f. correction (skip if G_sub <= 1)
					if (G_sub > 1) {
						V_slope = V_slope + sign * (G_sub / (G_sub - 1)) :* bread * meat_sub * bread
					}

					// Store cluster counts for single-variable subsets
					if (subset_size == 1) {
						for (j = 1; j <= n_clust; j++) {
							if (mask == 2^(j-1)) {
								st_numscalar("ivppmlhdfe_G" + strofreal(j), G_sub)
								if (G_sub < N_clust_min) N_clust_min = G_sub
							}
						}
					}
				}

				// Symmetrize
				_makesymmetric(V_slope)
				// PSD fix for multi-way clustering (matching reghdfe)
				if (n_clust > 1) {
					real scalar min_ev
					min_ev = min(symeigenvalues(V_slope))
					if (min_ev < 0) V_slope = V_slope - min_ev * I(K)
				}
			}
		}
	}
	else {
		V_slope = ivppmlhdfe_robustvce(Xhat, X_dm, irls_w, w_meat, resid, K, N_eff)
		{
			real matrix bread_r
			bread_r = invsym(cross(Xhat, irls_w, X_dm))
			actual_rank = K - diag0cnt(bread_r)
		}
	}

	// ---- Back-transform standardized coefficients and VCE ----
	// b_cons is invariant to standardization (mean(X_std)*b_std = mean(X)*b_raw)
	if (do_standardize & K > 0) {
		b = b :/ stdev_x'
		V_slope = V_slope :/ (stdev_x' * stdev_x)
		b_full = b \ b_cons
	}

	// Expand V to include _cons (VCE for _cons is zero: not estimable with absorbed FE)
	// Following ivreg2 partial() convention: partialled-out constant has no VCE
	V = J(K_total, K_total, 0)
	V[1..K, 1..K] = V_slope

	// ---- Log pseudo-likelihood ----
	// ll = sum(w * (y * eta - mu - lngamma(y + 1)))
	ll = quadsum(w_user :* (y :* eta - mu - lngamma(y :+ 1)))

	// Null model: intercept-only (offset-aware if offset present)
	// With offset: mu_0i = exp(c + o_i), c = log(sum(w*y) / sum(w*exp(o)))
	// Without offset: mu_0 = mean(y, w)
	{
		real scalar ll_0_c
		real colvector mu_0
		if (quadsum(abs(offset)) > 0) {
			ll_0_c = log(quadsum(w_user :* y) / quadsum(w_user :* exp(offset)))
			mu_0 = exp(ll_0_c :+ offset)
			ll_0 = quadsum(w_user :* (y :* (ll_0_c :+ offset) - mu_0 - lngamma(y :+ 1)))
		}
		else {
			ll_0_mu = mean(y, w_user)
			ll_0 = quadsum(w_user :* (y :* log(ll_0_mu) :- ll_0_mu :- lngamma(y :+ 1)))
		}
	}

	// ---- d() variable: save FE sum ----
	if (args() >= 17) {
		if (d_var_s != "") {
			real colvector d_vals
			d_vals = eta - offset
			if (K > 0) {
				d_vals = d_vals - X * b
			}
			d_vals = d_vals :- b_cons
			// Center to weighted mean zero (matching ppmlhdfe convention)
			d_vals = d_vals :- mean(d_vals, irls_w)
			// Use HDFE.sample (Stata obs indices) NOT touse_s — after a
			// separation trim the touse mask has more 1s than d_vals has
			// rows, causing st_store to crash with a conformability error.
			st_store(HDFE.sample, d_var_s, d_vals)
		}
	}

	// ---- Prepare HDFE for post() ----
	// Set fields that HDFE.post() reads (following ppmlhdfe approach)
	HDFE.depvar = depvar_s
	HDFE.indepvars = K > 0 ? (K_exog > 0 & K_endo > 0 ? (exog_vars, endog_vars) : (K_endo > 0 ? endog_vars : exog_vars)) : J(1, 0, "")
	HDFE.fullindepvars = invtokens(HDFE.indepvars)
	HDFE.title = "IV-PPML"
	HDFE.cmdline = "ivppmlhdfe " + st_local("0")
	HDFE.vcetype = vcetype_s
	HDFE.clustervars = has_cluster ? clust_varnames : J(1, 0, "")
	HDFE.base_clustervars = has_cluster ? clust_varnames : J(1, 0, "")
	HDFE.num_clusters = n_clust
	HDFE.compute_constant = 0
	HDFE.report_constant = 0
	HDFE.df_m = actual_rank
	HDFE.ll = ll
	HDFE.ll_0 = ll_0
	if (has_cluster) {
		HDFE.N_clust_list = J(1, n_clust, .)
		real scalar _j, _gmin
		_gmin = .
		for (_j = 1; _j <= n_clust; _j++) {
			HDFE.N_clust_list[_j] = st_numscalar("ivppmlhdfe_G" + strofreal(_j))
			if (HDFE.N_clust_list[_j] < _gmin) _gmin = HDFE.N_clust_list[_j]
		}
		HDFE.N_clust = _gmin
	}
	HDFE.estimate_dof()

	// Add _cons to indepvars for proper column naming
	HDFE.indepvars = HDFE.indepvars, "_cons"
	HDFE.fullindepvars = HDFE.fullindepvars + " _cons"
	HDFE.not_basevar = J(1, 0, .)
	{
		real scalar _kk
		for (_kk = 1; _kk <= K + 1; _kk++) {
			HDFE.not_basevar = HDFE.not_basevar, 1
		}
	}

	// Store HDFE as external so Stata-side code can call HDFE.post()
	{
		external class FixedEffects scalar ivppmlhdfe_HDFE
		ivppmlhdfe_HDFE = HDFE
	}

	// ---- Post to Stata ----
	st_matrix(bname, b_full')
	st_matrix(Vname, V)

	st_numscalar("ivppmlhdfe_converged", converged)
	st_numscalar("ivppmlhdfe_iterations", iter)
	st_numscalar("ivppmlhdfe_deviance", deviance)
	st_numscalar("ivppmlhdfe_N", N_eff)
	st_numscalar("ivppmlhdfe_N_rows", N)
	st_numscalar("ivppmlhdfe_ll", ll)
	st_numscalar("ivppmlhdfe_ll0", ll_0)
	st_numscalar("ivppmlhdfe_num_singletons", num_singletons)
	st_numscalar("ivppmlhdfe_actual_rank", actual_rank)
	st_numscalar("ivppmlhdfe_num_sep_fe", num_sep_fe)
	st_numscalar("ivppmlhdfe_num_sep_advanced", num_sep_advanced)
	st_numscalar("ivppmlhdfe_num_sep_mu", num_sep_mu)
}

end
