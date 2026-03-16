{smcl}
{* *! version 0.6.0  15mar2026}{...}
{vieweralsosee "ppmlhdfe" "help ppmlhdfe"}{...}
{vieweralsosee "reghdfe" "help reghdfe"}{...}
{vieweralsosee "ivreghdfe" "help ivreghdfe"}{...}
{viewerjumpto "Syntax" "ivppmlhdfe##syntax"}{...}
{viewerjumpto "Description" "ivppmlhdfe##description"}{...}
{viewerjumpto "Options" "ivppmlhdfe##options"}{...}
{viewerjumpto "Examples" "ivppmlhdfe##examples"}{...}
{viewerjumpto "Stored results" "ivppmlhdfe##results"}{...}
{viewerjumpto "References" "ivppmlhdfe##references"}{...}
{viewerjumpto "Author" "ivppmlhdfe##author"}{...}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col:{cmd:ivppmlhdfe} {hline 2}}IV Poisson pseudo-maximum likelihood estimation with high-dimensional fixed effects{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:ivppmlhdfe}
{depvar}
[{it:exogvars}]
{cmd:(}{it:endogvars} {cmd:=} {it:excluded_instruments}{cmd:)}
{ifin}
{weight}{cmd:,}
{cmdab:a:bsorb(}{it:absvars}{cmd:)}
[{it:options}]

{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{p2coldent:* {cmdab:a:bsorb(}{it:absvars}{cmd:)}}categorical variables to absorb as fixed effects{p_end}
{synopt:{cmdab:vce(}{it:vcetype}{cmd:)}}{it:vcetype} may be {cmd:robust} (default) or {cmd:cluster} {it:clustvar}{p_end}
{synopt:{cmdab:cl:uster(}{it:clustvar}{cmd:)}}synonym for {cmd:vce(cluster} {it:clustvar}{cmd:)}{p_end}
{synopt:{cmdab:tol:erance(}{it:#}{cmd:)}}convergence tolerance; default is {cmd:1e-8}{p_end}
{synopt:{cmdab:maxit:erations(}{it:#}{cmd:)}}maximum number of IRLS iterations; default is {cmd:1000}{p_end}
{synopt:{cmdab:v:erbose(}{it:#}{cmd:)}}verbosity level 0/1/2/3; default is {cmd:0}{p_end}
{synopt:{cmd:nolog}}suppress iteration log{p_end}
{synopt:{cmdab:ef:orm}}display exponentiated coefficients{p_end}
{synopt:{cmd:irr}}synonym for {cmd:eform}{p_end}
{synopt:{cmdab:bias:correction(}{it:class id1 id2}{cmd:)}}analytical bias correction; see {help ivppmlhdfe##biascorrection:below}{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}* {cmd:absorb()} is required.{p_end}
{p 4 6 2}{cmd:pweight}s and {cmd:fweight}s are allowed.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:ivppmlhdfe} estimates an instrumental-variables Poisson pseudo-maximum
likelihood (IV-PPML) model with high-dimensional fixed effects.  The estimator
uses iteratively reweighted GMM following
{help ivppmlhdfe##references:Mullahy (1997)}, with fixed effects concentrated
out via the PPML first-order condition at each iteration using
{help reghdfe}.

{pstd}
The command handles models with multiple sets of high-dimensional fixed
effects, including the two-way (exporter-year, importer-year) and three-way
(exporter-year, importer-year, pair) structures common in gravity models
of international trade.

{pstd}
Unlike {cmd:ivpoisson} (which supports only simple factor-variable fixed effects),
{cmd:ivppmlhdfe} can absorb hundreds of thousands of fixed effects efficiently
using the iterative demeaning algorithm of {cmd:reghdfe}.


{marker options}{...}
{title:Options}

{phang}
{cmd:absorb(}{it:absvars}{cmd:)} specifies the categorical variables to be
absorbed as fixed effects.  Multiple variables are separated by spaces.
Interactions are specified using {cmd:#} notation as in {cmd:reghdfe}.
This option is required.

{phang}
{cmd:vce(robust)} (the default) computes the Eicker-Huber-White
heteroskedasticity-robust sandwich variance estimator.
{cmd:vce(cluster} {it:clustvar}{cmd:)} computes the cluster-robust
variance estimator.

{phang}
{cmd:tolerance(}{it:#}{cmd:)} sets the convergence criterion for the
IRLS algorithm.  Default is {cmd:1e-8}.  Convergence is declared when the
maximum absolute change in coefficients falls below this threshold.

{phang}
{cmd:maxiterations(}{it:#}{cmd:)} sets the maximum number of IRLS
iterations.  Default is {cmd:1000}.

{phang}
{cmd:verbose(}{it:#}{cmd:)} controls the amount of output displayed
during estimation.  0 = quiet (default), 1 = iteration summary,
2 = detailed iteration info, 3 = full debug output.

{phang}
{cmd:eform} (or equivalently {cmd:irr}) displays exponentiated
coefficients, interpretable as incidence-rate ratios.

{marker biascorrection}{...}
{phang}
{cmd:biascorrection(}{it:class id1 id2}{cmd:)} requests an analytical
incidental-parameter bias correction for the point estimates.
{it:class} specifies the fixed-effect structure:

{phang2}
{cmd:biascorrection(a} {it:id_var time_var}{cmd:)} — Class A model with
individual ({it:id_var}) and time ({it:time_var}) fixed effects.
Corrects for O(1/T) + O(1/N) bias.

{phang2}
{cmd:biascorrection(b} {it:exp_time_var imp_time_var}{cmd:)} — Class B
gravity model with exporter-time and importer-time fixed effects.
Corrects for O(1/N) bias.

{phang2}
Requires {cmd:vce(cluster} {it:clustvar}{cmd:)}.


{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse trade}{p_end}

{pstd}Basic IV-PPML with time fixed effects{p_end}
{phang2}{cmd:. ivppmlhdfe trade exog_var (endog_var = instrument), absorb(year) vce(robust)}{p_end}

{pstd}Two-way gravity model (exporter-year, importer-year FE){p_end}
{phang2}{cmd:. ivppmlhdfe trade (policy = instrument), absorb(exp_year imp_year) vce(cluster pair_id)}{p_end}

{pstd}Three-way gravity model with pair FE{p_end}
{phang2}{cmd:. ivppmlhdfe trade (policy = instrument), absorb(exp_year imp_year pair_id) vce(cluster pair_id)}{p_end}

{pstd}With bias correction (Class A){p_end}
{phang2}{cmd:. ivppmlhdfe y (x = z), absorb(id year) vce(cluster id) biascorrection(a id year)}{p_end}

{pstd}Display incidence-rate ratios{p_end}
{phang2}{cmd:. ivppmlhdfe trade (policy = instrument), absorb(exp_year imp_year) vce(cluster pair_id) irr}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:ivppmlhdfe} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(iterations)}}number of IRLS iterations{p_end}
{synopt:{cmd:e(converged)}}1 if converged, 0 otherwise{p_end}
{synopt:{cmd:e(ll)}}log pseudolikelihood (Poisson){p_end}
{synopt:{cmd:e(deviance)}}deviance{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:ivppmlhdfe}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(exog)}}exogenous regressors{p_end}
{synopt:{cmd:e(endog)}}endogenous regressors{p_end}
{synopt:{cmd:e(instruments)}}excluded instruments{p_end}
{synopt:{cmd:e(absorb)}}absorbed fixed-effect variables{p_end}
{synopt:{cmd:e(vcetype)}}variance estimator type{p_end}
{synopt:{cmd:e(clustvar)}}cluster variable (if applicable){p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}


{marker references}{...}
{title:References}

{phang}
Mullahy, J. 1997.
Instrumental-variable estimation of count data models: Applications to
models of cigarette smoking behavior.
{it:Review of Economics and Statistics} 79(4): 586-593.

{phang}
Correia, S., P. Guimaraes, and T. Zylkin. 2020.
Fast Poisson estimation with high-dimensional fixed effects.
{it:Stata Journal} 20(1): 95-115.

{phang}
Weidner, M. and T. Zylkin. 2021.
Bias and consistency in three-way gravity models.
{it:Journal of International Economics} 132: 103513.

{phang}
Fernandez-Val, I. and M. Weidner. 2016.
Individual and time effects in nonlinear panel models with large N, T.
{it:Journal of Econometrics} 192(1): 291-312.

{phang}
Kauermann, G. and R. J. Carroll. 2001.
A note on the efficiency of sandwich covariance matrix estimation.
{it:Journal of the American Statistical Association} 96(456): 1387-1396.


{marker author}{...}
{title:Author}

{pstd}Ohyun Kwon{p_end}
