*! ivppmlhdfe predict program
*! Follows ppmlhdfe_p.ado pattern

program ivppmlhdfe_p
	version 14.0

	* Intercept -scores- for margins compatibility
	cap syntax anything [if] [in], SCores
	if !c(rc) {
		_score_spec `anything', score
		local 0 `s(varlist)' `if' `in', scores
	}

	syntax newvarname [if] [in] [, ///
		XB XBD D Mu N Eta STDP ///
		SCores RESiduals RESPonse PEARson Anscombe WORKing DEViance]

	* Count how many options specified
	local opts `xb' `xbd' `d' `mu' `n' `eta' `stdp' `scores' `residuals' `response' `pearson' `anscombe' `working' `deviance'
	opts_exclusive "`opts'"
	local opt `opts'

	* Default: mu
	if "`opt'" == "" {
		di as txt "(option mu assumed; predicted mean of depvar)"
		local opt "mu"
	}

	* Alias: n -> mu, eta -> xbd
	if "`opt'" == "n" local opt "mu"
	if "`opt'" == "eta" local opt "xbd"

	* stdp delegates to _predict on the slope-block VCE.
	* Note: e(V) for _cons is zero by design (partialled-out, ivreg2 convention),
	* so stdp reflects only the slope-block contribution at FE centers.
	if "`opt'" == "stdp" {
		_predict double `varlist' `if' `in', stdp
		label var `varlist' "S.E. of linear prediction"
		exit
	}

	* Check if d() is needed (all options except xb/stdp)
	if !inlist("`opt'", "xb") {
		if "`e(absvars)'" != "_cons" {
			_assert `"`e(d)'"' != "", ///
				msg("predict `opt' requires the {bf:d()} option during estimation")
			confirm double var `e(d)', exact
		}
	}

	* Compute xb = X*b + _cons
	tempvar xbvar
	PredictXB double `xbvar' `if' `in'

	* For xb, we are done
	if "`opt'" == "xb" {
		gen double `varlist' = `xbvar' `if' `in'
		label var `varlist' "Linear prediction (xb)"
		exit
	}

	* Compute full linear predictor: eta = xb + d
	* Note: _predict (via PredictXB) already includes e(offset) when set,
	* and d = eta - offset - X*b - cons, so xb + d = X*b + cons + offset + FE = eta.
	* Do NOT add e(offset) again here (would double-count).
	tempvar etavar
	if "`e(absvars)'" == "_cons" {
		* Constant-only model: no FE
		gen double `etavar' = `xbvar' `if' `in'
	}
	else {
		gen double `etavar' = `xbvar' + `e(d)' `if' `in'
	}

	* d (fixed effects only)
	if "`opt'" == "d" {
		if "`e(absvars)'" == "_cons" {
			gen double `varlist' = 0 `if' `in'
		}
		else {
			gen double `varlist' = `e(d)' `if' `in'
		}
		label var `varlist' "Sum of fixed effects"
		exit
	}

	* xbd / eta (full linear predictor)
	if "`opt'" == "xbd" {
		gen double `varlist' = `etavar' `if' `in'
		label var `varlist' "Linear prediction (xb + d)"
		exit
	}

	* mu (predicted count)
	if "`opt'" == "mu" {
		gen double `varlist' = exp(`etavar') `if' `in'
		label var `varlist' "Predicted mean (mu)"
		exit
	}

	* scores / response / residuals: y - mu
	if inlist("`opt'", "scores", "residuals", "response") {
		tempvar muvar
		gen double `muvar' = exp(`etavar') `if' `in'
		gen double `varlist' = `e(depvar)' - `muvar' `if' `in'
		if "`opt'" == "scores" label var `varlist' "Score (y - mu)"
		else label var `varlist' "Response residual (y - mu)"
		exit
	}

	* pearson: (y - mu) / sqrt(mu)
	if "`opt'" == "pearson" {
		tempvar muvar
		gen double `muvar' = exp(`etavar') `if' `in'
		gen double `varlist' = (`e(depvar)' - `muvar') / sqrt(`muvar') `if' `in'
		label var `varlist' "Pearson residual"
		exit
	}

	* anscombe: 1.5 * (y^(2/3) - mu^(2/3)) / mu^(1/6)  (matching ppmlhdfe)
	if "`opt'" == "anscombe" {
		tempvar muvar
		gen double `muvar' = exp(`etavar') `if' `in'
		gen double `varlist' = 1.5 * (`e(depvar)'^(2/3) - `muvar'^(2/3)) / `muvar'^(1/6) `if' `in'
		label var `varlist' "Anscombe residual"
		exit
	}

	* working: (y - mu) / mu
	if "`opt'" == "working" {
		tempvar muvar
		gen double `muvar' = exp(`etavar') `if' `in'
		gen double `varlist' = (`e(depvar)' - `muvar') / `muvar' `if' `in'
		label var `varlist' "Working residual"
		exit
	}

	* deviance contribution: 2 * (mu - y + y*ln(y/mu))  (matching ppmlhdfe)
	if "`opt'" == "deviance" {
		tempvar muvar
		gen double `muvar' = exp(`etavar') `if' `in'
		gen double `varlist' = 2 * cond(`e(depvar)'>0, ///
			`muvar' - `e(depvar)' + `e(depvar)'*ln(`e(depvar)'/`muvar'), ///
			`muvar') `if' `in'
		label var `varlist' "Deviance contribution"
		exit
	}

	di as err "option `opt' not recognized"
	exit 198
end


* Helper: compute xb = X*b + _cons
* Handles the edge case where no regressors exist
program PredictXB
	syntax newvarname [if] [in] [, *]
	capture _predict `typlist' `varlist' `if' `in', `options'
	if c(rc) {
		* No regressors: xb = _cons
		gen double `varlist' = _b[_cons] `if' `in'
	}
end
