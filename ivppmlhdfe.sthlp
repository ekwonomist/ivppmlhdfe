{smcl}
{* *! version 0.9.4  23apr2026}{...}
{vieweralsosee "ppmlhdfe" "help ppmlhdfe"}{...}
{vieweralsosee "reghdfe" "help reghdfe"}{...}
{vieweralsosee "ivreghdfe" "help ivreghdfe"}{...}
{viewerjumpto "Syntax" "ivppmlhdfe##syntax"}{...}
{viewerjumpto "Description" "ivppmlhdfe##description"}{...}
{viewerjumpto "Options" "ivppmlhdfe##options"}{...}
{viewerjumpto "Examples" "ivppmlhdfe##examples"}{...}
{viewerjumpto "Predict" "ivppmlhdfe##predict"}{...}
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
{cmd:(}{it:endogvar} {cmd:=} {it:excluded_instrument}{cmd:)}
{ifin}
{weight}{cmd:,}
{cmdab:a:bsorb(}{it:absvars}{cmd:)}
[{it:options}]

{synoptset 36 tabbed}{...}
{synopthdr}
{synoptline}
{p2coldent:* {cmdab:a:bsorb(}{it:absvars}{cmd:)}}categorical variables to absorb as fixed effects{p_end}
{synopt:{cmdab:vce(}{it:vcetype}{cmd:)}}{it:vcetype} may be {cmd:robust} (default) or {cmd:cluster} {it:c1} [{it:c2} [{it:c3}]]{p_end}
{synopt:{cmdab:exp:osure(}{it:varname}{cmd:)}}exposure variable; log({it:varname}) used as offset{p_end}
{synopt:{cmdab:off:set(}{it:varname}{cmd:)}}offset variable with coefficient fixed at 1{p_end}
{synopt:{cmdab:d(}{it:newvarname}{cmd:)}}save sum of absorbed fixed effects{p_end}
{synopt:{cmd:d2}}equivalent to {cmd:d(_ivppmlhdfe_d)} (auto-named){p_end}
{synopt:{cmdab:sep:aration(}{it:method}{cmd:)}}separation detection: {cmd:fe simplex relu} (default), {cmd:all} (adds {cmd:mu}), or {cmd:none}{p_end}
{synopt:{cmdab:stand:ardize}}standardize X and Z columns to unit variance for numerical stability{p_end}
{synopt:{cmdab:guess(}{it:method}{cmd:)}}initial mu guess: {cmd:simple} (default) or {cmd:mean}{p_end}
{synopt:{cmdab:tol:erance(}{it:#}{cmd:)}}outer convergence tolerance; default is {cmd:1e-8}{p_end}
{synopt:{cmdab:itol:erance(}{it:#}{cmd:)}}inner FE-solver target tolerance; default is {cmd:max(1e-12, 0.1*tol)}{p_end}
{synopt:{cmdab:maxit:erations(}{it:#}{cmd:)}}maximum IRLS iterations; default is {cmd:10000}{p_end}
{synopt:{cmdab:keepsin:gletons}}keep singleton observations (default drops them){p_end}
{synopt:{cmdab:tagsep(}{it:newvarname}{cmd:)}}separation-detection-only mode: tag separated obs and exit{p_end}
{synopt:{cmdab:zvar:name(}{it:newvarname}{cmd:)}}save ReLU separation certificate variable{p_end}
{synopt:{cmdab:v:erbose(}{it:#}{cmd:)}}verbosity level {cmd:-1}/0/1/2; default is {cmd:0}{p_end}
{synopt:{cmd:nolog}}suppress iteration log (results table still shown){p_end}
{synopt:{cmdab:ef:orm}}display exponentiated coefficients{p_end}
{synopt:{cmd:irr}}synonym for {cmd:eform}{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}* {cmd:absorb()} is required (or specify {cmd:noabsorb} for no fixed effects).{p_end}
{p 4 6 2}{cmd:pweight}s and {cmd:fweight}s are allowed.{p_end}
{p 4 6 2}Both {bf:just-identified} and {bf:overidentified} models are supported. Multiple endogenous regressors are allowed.{p_end}
{p 4 6 2}Unrecognised options (e.g. {cmd:dof()}, {cmd:pool()}, {cmd:accel()}) are forwarded to {cmd:reghdfe} via the {cmd:absorb()} string.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:ivppmlhdfe} estimates an instrumental-variable Poisson pseudo-maximum
likelihood (IV-PPML) model with high-dimensional fixed effects. The estimator
targets the additive moment condition E[q(y - mu)] = 0 with
q = (x', z')' (exogenous regressors stacked with excluded instruments),
together with the per-group fixed-effect score sum_{g in r}(y_g - mu_g) = 0,
following {help ivppmlhdfe##references:Windmeijer and Santos Silva (1997)}.
It is solved via iteratively reweighted 2SLS
({help ivppmlhdfe##references:Correia, Guimaraes, and Zylkin 2020}), with
fixed effects concentrated out at each iteration using {help reghdfe}.

{pstd}
The command handles models with multiple sets of high-dimensional fixed
effects, including the two-way (exporter-year, importer-year) and three-way
(exporter-year, importer-year, pair) structures common in gravity models
of international trade.

{pstd}
For {bf:just-identified} models, the 2SLS first-order condition coincides
with the moment system above. For {bf:overidentified} models, the 2SLS
imposes a projected condition (analogous to {help ivreg2}).

{pstd}
{bf:Incidental parameter problem (IPP) and the SPJ remedy.} The Poisson
'score cancellation' that renders standard PPML immune to leading-order
incidental parameter bias does not carry over to the IV-PPML case, so
estimates here can be biased even when standard PPML is not. For the bias
orders by FE structure (individual + time, source-time + dest-time,
source-time + dest-time + pair), the split-panel jackknife (SPJ)
formulas, the recommended bootstrap aggregator, and Monte Carlo evidence,
see Kwon, Larch, Yoon, and Yotov (2026). Ready-to-run SPJ + bootstrap
templates for each FE class are shipped under {bf:data/} on the GitHub
repository.

{pstd}
v0.9.4 automatically detects and drops regressors and instruments that are
collinear with each other or fully absorbed by the fixed effects.  Dropped
variables appear as omitted in the coefficient table.

{pstd}
The command is designed to feel natural to users of {cmd:ppmlhdfe}.
All standard post-estimation commands are supported: {cmd:test}, {cmd:lincom},
{cmd:nlcom}, {cmd:predict}, {cmd:margins}, {cmd:estat}, {cmd:regsave},
{cmd:esttab}, {cmd:outreg2}, {cmd:coefplot}.


{marker options}{...}
{title:Options}

{phang}
{cmd:absorb(}{it:absvars}{cmd:)} specifies the categorical variables to be
absorbed as fixed effects.  Multiple variables are separated by spaces.
Interactions are specified using {cmd:#} notation as in {cmd:reghdfe}.
Factor variable syntax (e.g., {cmd:i.id#i.year}) is supported.
This option is required.

{phang}
{cmd:vce(robust)} (the default) computes the Eicker-Huber-White
heteroskedasticity-robust sandwich variance estimator with N/(N-K)
small-sample correction.
{cmd:vce(cluster} {it:c1} [{it:c2} [{it:c3}]]{cmd:)} computes the
cluster-robust variance estimator.  Multi-way clustering uses the
Cameron-Gelbach-Miller (2011) inclusion-exclusion formula with G/(G-1)
correction per cluster dimension.

{phang}
{cmd:exposure(}{it:varname}{cmd:)} specifies an exposure variable.
The log of {it:varname} is included as an offset with coefficient
constrained to 1.  Mutually exclusive with {cmd:offset()}.

{phang}
{cmd:offset(}{it:varname}{cmd:)} specifies an offset variable included
in the linear predictor with coefficient constrained to 1.
Mutually exclusive with {cmd:exposure()}.

{phang}
{cmd:d(}{it:newvarname}{cmd:)} saves the sum of absorbed fixed effects
as a new variable.  Required for {cmd:predict mu} and most other
predict statistics (except {cmd:predict xb}).

{phang}
{cmd:d2} is shorthand for {cmd:d(_ivppmlhdfe_d)}; ivppmlhdfe will create
the variable {cmd:_ivppmlhdfe_d} (dropping any pre-existing one).

{phang}
{cmd:separation(}{it:method list}{cmd:)} controls the techniques used to
detect and drop separated observations. Valid techniques are
{cmd:fe simplex relu mu}; aliases {cmd:default} (or omitted) =
{cmd:fe simplex relu}, {cmd:all} = {cmd:fe simplex relu mu}, and
{cmd:none} disables detection. {cmd:fe} drops fixed-effect groups where
all dependent-variable observations are zero (via the {cmd:reghdfe}
iweight trick, which correctly handles interaction FEs such as
{cmd:i.exp#i.imp}). {cmd:simplex} and {cmd:relu} use {cmd:ppmlhdfe}'s
linear-programming detection. {cmd:mu} flags observations whose linear
predictor diverges to negative infinity at run-time. Singletons are
dropped by default via {cmd:reghdfe}.

{phang}
{cmd:standardize} divides each X and Z column by its weighted standard
deviation before estimation, which can help numerical stability with
poorly-scaled regressors. Coefficients and the variance matrix are
back-transformed after convergence and are identical to the
unstandardized fit (up to floating-point precision).

{phang}
{cmd:guess(}{it:method}{cmd:)} controls the initial value of mu used to
start the IRLS loop. {cmd:simple} (the default) uses
mu_0 = 0.5*(y + mean(y)); {cmd:mean} uses mu_0 = mean(y) for every
observation.

{phang}
{cmd:tolerance(}{it:#}{cmd:)} sets the outer convergence criterion for
the IRLS algorithm.  Default is {cmd:1e-8}.  Convergence is declared
when the relative change in deviance falls below this threshold.

{phang}
{cmd:itolerance(}{it:#}{cmd:)} sets the inner FE-solver target tolerance
used by the adaptive scheme. Default is {cmd:max(1e-12, 0.1 * tolerance)}.
Tighten this if you suspect under-converged demeaning at the inner step.

{phang}
{cmd:maxiterations(}{it:#}{cmd:)} sets the maximum number of IRLS
iterations.  Default is {cmd:10000}.

{phang}
{cmd:keepsingletons} prevents singleton observations from being dropped
by {cmd:reghdfe}. Use with caution; cluster-robust standard errors are
not generally valid in the presence of singletons.

{phang}
{cmd:tagsep(}{it:newvarname}{cmd:)} runs separation detection only (no
IRLS) and saves a 0/1 indicator for separated observations to
{it:newvarname}, then exits without estimating coefficients. Useful
for diagnosing which observations a downstream estimator would drop.

{phang}
{cmd:zvarname(}{it:newvarname}{cmd:)} saves the ReLU separation
certificate variable to {it:newvarname}; only meaningful when
{cmd:relu} is in the {cmd:separation()} list.

{phang}
{cmd:verbose(}{it:#}{cmd:)} controls the amount of output displayed
during estimation.  -1 = silent, 0 = summary only (default),
1 = iteration log, 2 = detailed output.

{phang}
{cmd:nolog} suppresses the IRLS iteration log but still displays the
results table.

{phang}
{cmd:eform} (or equivalently {cmd:irr}) displays exponentiated
coefficients, interpretable as incidence-rate ratios.


{marker examples}{...}
{title:Examples}

{pstd}The examples below use the three datasets shipped with the package in
the {cmd:data/} folder of the GitHub repository: {cmd:ivppmlhdfe_ClassA.dta}
(individual + time FE), {cmd:ivppmlhdfe_ClassB.dta} (two-way gravity FE), and
{cmd:ivppmlhdfe_ClassC.dta} (three-way gravity FE). Each has dependent
variable {cmd:y}, exogenous regressor {cmd:x2}, endogenous regressor {cmd:x1},
and excluded instrument {cmd:z}.{p_end}

{pstd}Class A: individual + time FE{p_end}
{phang2}{cmd:. use "ivppmlhdfe_ClassA.dta", clear}{p_end}
{phang2}{cmd:. ivppmlhdfe y x2 (x1 = z), absorb(id year) vce(robust)}{p_end}

{pstd}Class B: two-way gravity FE (exporter-year, importer-year){p_end}
{phang2}{cmd:. use "ivppmlhdfe_ClassB.dta", clear}{p_end}
{phang2}{cmd:. ivppmlhdfe y x2 (x1 = z), absorb(exp#year imp#year) vce(cluster pair)}{p_end}

{pstd}Class C: three-way gravity FE (pair, exporter-year, importer-year){p_end}
{phang2}{cmd:. use "ivppmlhdfe_ClassC.dta", clear}{p_end}
{phang2}{cmd:. ivppmlhdfe y x2 (x1 = z), absorb(exp#imp exp#year imp#year) vce(cluster pair)}{p_end}

{pstd}Multi-way clustering (exporter and importer){p_end}
{phang2}{cmd:. use "ivppmlhdfe_ClassB.dta", clear}{p_end}
{phang2}{cmd:. ivppmlhdfe y x2 (x1 = z), absorb(exp#year imp#year) vce(cluster exp imp)}{p_end}

{pstd}With {cmd:d()} for prediction{p_end}
{phang2}{cmd:. ivppmlhdfe y x2 (x1 = z), absorb(exp#year imp#year) d(fe_sum)}{p_end}
{phang2}{cmd:. predict muhat, mu}{p_end}
{phang2}{cmd:. predict resid, residuals}{p_end}

{pstd}Display incidence-rate ratios{p_end}
{phang2}{cmd:. ivppmlhdfe y x2 (x1 = z), absorb(exp#year imp#year) vce(cluster pair) irr}{p_end}

{pstd}Margins{p_end}
{phang2}{cmd:. ivppmlhdfe y x2 (x1 = z), absorb(id year) d(fe_sum)}{p_end}
{phang2}{cmd:. margins, dydx(x2)}{p_end}

{pstd}Standardize regressors for numerical stability{p_end}
{phang2}{cmd:. ivppmlhdfe y x2 (x1 = z), absorb(id year) standardize}{p_end}

{pstd}Aggressive separation detection (FE + simplex + ReLU + mu){p_end}
{phang2}{cmd:. ivppmlhdfe y x2 (x1 = z), absorb(id year) separation(all)}{p_end}

{pstd}Tag separated observations without estimating{p_end}
{phang2}{cmd:. ivppmlhdfe y x2 (x1 = z), absorb(id year) separation(all) tagsep(sep_obs)}{p_end}


{marker predict}{...}
{title:Predict}

{pstd}
After estimation, {cmd:predict} supports the following statistics:

{synoptset 20 tabbed}{...}
{synopthdr:statistic}
{synoptline}
{synopt:{cmd:mu}}predicted mean exp(eta); the default{p_end}
{synopt:{cmd:n}}synonym for {cmd:mu}{p_end}
{synopt:{cmd:xb}}linear prediction X*b + _cons (includes offset if present){p_end}
{synopt:{cmd:xbd} or {cmd:eta}}full linear predictor xb + d{p_end}
{synopt:{cmd:d}}sum of absorbed fixed effects{p_end}
{synopt:{cmd:scores}}score residual y - mu{p_end}
{synopt:{cmd:residuals}}response residual y - mu{p_end}
{synopt:{cmd:pearson}}Pearson residual (y - mu) / sqrt(mu){p_end}
{synopt:{cmd:working}}working residual (y - mu) / mu{p_end}
{synopt:{cmd:deviance}}deviance residual{p_end}
{synopt:{cmd:stdp}}standard error of the linear prediction (slope block only; _cons has no SE by design){p_end}
{synopt:{cmd:anscombe}}Anscombe residual: 1.5*(y^(2/3) - mu^(2/3)) / mu^(1/6){p_end}
{synoptline}

{pstd}
All statistics except {cmd:xb} require the {cmd:d()} option during estimation.
{cmd:margins} works automatically via the {cmd:scores} option.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:ivppmlhdfe} stores the following in {cmd:e()}:

{synoptset 28 tabbed}{...}
{p2col 5 28 32 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_full)}}observations before singleton/separation drops{p_end}
{synopt:{cmd:e(num_separated)}}total separated observations dropped{p_end}
{synopt:{cmd:e(num_sep_fe)}}separated by {cmd:fe} method{p_end}
{synopt:{cmd:e(num_sep_advanced)}}separated by {cmd:simplex}/{cmd:relu} methods{p_end}
{synopt:{cmd:e(num_sep_mu)}}separated by {cmd:mu} method (in-loop){p_end}
{synopt:{cmd:e(num_singletons)}}number of singletons dropped{p_end}
{synopt:{cmd:e(iterations)}}number of IRLS iterations{p_end}
{synopt:{cmd:e(converged)}}1 if converged, 0 otherwise{p_end}
{synopt:{cmd:e(deviance)}}deviance at convergence{p_end}
{synopt:{cmd:e(ll)}}log pseudo-likelihood (descriptive){p_end}
{synopt:{cmd:e(ll_0)}}null-model log pseudo-likelihood{p_end}
{synopt:{cmd:e(r2_p)}}pseudo R-squared (1 - ll/ll_0){p_end}
{synopt:{cmd:e(chi2)}}Wald chi-squared statistic{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom (actual rank){p_end}
{synopt:{cmd:e(rank)}}rank of e(V){p_end}
{synopt:{cmd:e(N_clustervars)}}number of cluster variables{p_end}
{synopt:{cmd:e(N_clust)}}minimum number of clusters{p_end}
{synopt:{cmd:e(N_clust1)}}number of clusters (1st variable){p_end}
{synopt:{cmd:e(N_clust2)}}number of clusters (2nd variable){p_end}

{synoptset 28 tabbed}{...}
{p2col 5 28 32 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:ivppmlhdfe}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(exogvars)}}exogenous regressors{p_end}
{synopt:{cmd:e(endogvars)}}endogenous regressors{p_end}
{synopt:{cmd:e(instruments)}}excluded instruments{p_end}
{synopt:{cmd:e(absorb)}}absorb() specification{p_end}
{synopt:{cmd:e(absvars)}}absorbed fixed-effect variables{p_end}
{synopt:{cmd:e(vcetype)}}variance estimator type{p_end}
{synopt:{cmd:e(clustvar)}}cluster variable(s){p_end}
{synopt:{cmd:e(title)}}title for estimation output{p_end}
{synopt:{cmd:e(chi2type)}}type of chi-squared test{p_end}
{synopt:{cmd:e(predict)}}predict program name{p_end}
{synopt:{cmd:e(marginsok)}}predict statistics allowed by margins{p_end}
{synopt:{cmd:e(properties)}}estimation properties{p_end}
{synopt:{cmd:e(separation)}}separation detection method{p_end}
{synopt:{cmd:e(offset)}}offset variable name{p_end}
{synopt:{cmd:e(d)}}d() variable name{p_end}

{synoptset 28 tabbed}{...}
{p2col 5 28 32 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector (slopes + _cons){p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix{p_end}

{synoptset 28 tabbed}{...}
{p2col 5 28 32 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}

{pstd}
{bf:Note:} The {cmd:_cons} entry in {cmd:e(V)} is zero.  The constant is a
derived quantity (average of absorbed fixed effects) with no estimable
variance, following {cmd:ivreg2}'s convention for partialled-out variables.

{pstd}
{bf:Note:} {cmd:e(ll)} and {cmd:e(ll_0)} are descriptive quantities evaluated
at the IV solution.  They are not the criterion being optimized, so
likelihood-ratio tests are not valid with endogenous regressors.

{pstd}
{bf:Error codes:}

{synoptset 28 tabbed}{...}
{p2col 5 28 32 2: Return codes}{p_end}
{synopt:{cmd:rc=430}}Non-convergence.  The IRLS iteration failed to reach
tolerance within {cmd:maxiterations()} iterations.  Try increasing
{cmd:maxiterations()}, or use {cmd:standardize} or {cmd:mu} separation
options.{p_end}

{pstd}
{bf:Known limitation:} {cmd:_cons} is reported without a standard error.
This is by design: the constant is recovered from the weighted mean of
the partialled-out linear predictor (following the {cmd:ivreg2} convention)
and does not have an analytical SE.  Use bootstrap for inference on
{cmd:_cons}.


{marker references}{...}
{title:References}

{phang}
Windmeijer, F. A. G., and J. M. C. Santos Silva. 1997.
Endogeneity in count data models: An application to demand for health care.
{it:Journal of Applied Econometrics} 12(3): 281-294.

{phang}
Correia, S., P. Guimaraes, and T. Zylkin. 2020.
Fast Poisson estimation with high-dimensional fixed effects.
{it:Stata Journal} 20(1): 95-115.

{phang}
Kwon, O., M. Larch, J. Yoon, and Y. V. Yotov. 2026.
Instrumental-Variable Poisson PML with High-Dimensional Fixed Effects.
Working paper. (For IPP bias orders, the SPJ + bootstrap remedy, and
Monte Carlo evidence.)

{phang}
Weidner, M. and T. Zylkin. 2021.
Bias and consistency in three-way gravity models.
{it:Journal of International Economics} 132: 103513.

{phang}
Fernandez-Val, I. and M. Weidner. 2016.
Individual and time effects in nonlinear panel models with large N, T.
{it:Journal of Econometrics} 192(1): 291-312.

{phang}
Cameron, A. C., J. B. Gelbach, and D. L. Miller. 2011.
Robust inference with multiway clustering.
{it:Journal of Business and Economic Statistics} 29(2): 238-249.

{phang}
Mullahy, J. 1997.
Instrumental-variable estimation of count data models: Applications to
models of cigarette smoking behavior.
{it:Review of Economics and Statistics} 79(4): 586-593.
(Related but distinct: a transformation estimator for multiplicative
unobserved heterogeneity, not the additive moment used here.)


{marker author}{...}
{title:Authors}

{pstd}Ohyun Kwon, Mario Larch, Jangsu Yoon, Yoto V. Yotov{p_end}

{pstd}Report issues to Ohyun Kwon at {browse "mailto:theekwonomist@gmail.com":theekwonomist@gmail.com} or open a {browse "https://github.com/ekwonomist/ivppmlhdfe/issues":GitHub issue}.{p_end}
