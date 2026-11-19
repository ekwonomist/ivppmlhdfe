\name{ivppmlhdfe}
\alias{ivppmlhdfe}
\title{IV-PPML with high-dimensional fixed effects}
\description{
Estimates a Poisson pseudo-maximum likelihood model with endogenous regressors, excluded
instruments and any number of (interacted) high-dimensional fixed effects. Without an IV part
the model is plain PPML and reproduces \code{ppmlhdfe}.
}
\usage{
ivppmlhdfe(fml, data, vcov = "robust", weights = NULL,
           weight_type = c("pweight", "fweight"), offset = NULL, exposure = NULL,
           accel = TRUE, ffirst = FALSE, tol = 1e-8, itol = NULL, maxit = 10000L,
           separation = c("fe", "relu"), keep_singletons = FALSE, verbose = 0L)
}
\arguments{
  \item{fml}{a formula \code{y ~ exog | fe | endog ~ instruments}. Parts are separated by
    \code{|}: the exogenous regressors, the fixed effects (terms joined by \code{+},
    interactions written \code{exp^year}), and the IV part. \code{y ~ exog | fe} is plain PPML,
    \code{y ~ exog | endog ~ instruments} is IV without fixed effects, \code{y ~ exog} is
    PPML with a constant only. Regressor terms are expanded with \code{model.matrix}
    (factors become dummies; the intercept is absorbed).}
  \item{data}{a data frame.}
  \item{vcov}{\code{"robust"} (heteroskedasticity-robust, the default), a one-sided formula of
    cluster variables such as \code{~pair} (one-way) or \code{~exp + imp} (multi-way, Cameron,
    Gelbach and Miller), or a character vector of cluster variable names. Interacted clusters
    are written \code{~exp^imp}.}
  \item{weights}{\code{NULL} or the name of a weight variable (or a numeric vector).
    Observations with zero weight are dropped; negative weights are an error.}
  \item{weight_type}{\code{"pweight"} (sampling weights) or \code{"fweight"} (frequency weights:
    the reported N is the sum of the weights and the robust meat scales as if rows were
    replicated).}
  \item{offset}{\code{NULL} or the name of an offset variable (coefficient fixed at 1).}
  \item{exposure}{\code{NULL} or the name of an exposure variable (offset = log(exposure)).}
  \item{accel}{\code{TRUE} (default) switches on the fast path of version 0.10.0: (w) warm
    starts of the fixed-effect partial-out across IRLS iterations, (t) the staged inner
    tolerance schedule of fereg, and (r) a refresh of the partialled columns at the converged
    weights for the standard errors and first-stage statistics. \code{FALSE} gives the
    reference path: cold-started partial-outs, the inner-tolerance rule of ppmlhdfe and the
    standard-error convention of version 0.9.4. A string of letters (for example \code{"wt"})
    switches components on one by one (diagnostics). The fixed-effect solver (exact
    elimination of the largest fixed effect, then Jacobi-preconditioned conjugate gradients)
    is used in all cases.}
  \item{ffirst}{logical: report first-stage statistics for each endogenous regressor and the
    identification tests (see Value).}
  \item{tol}{convergence tolerance of the IRLS loop (relative change in the deviance).}
  \item{itol}{target tolerance of the inner fixed-effect solver; default \code{max(1e-12, 0.1 * tol)}.}
  \item{maxit}{maximum number of IRLS iterations.}
  \item{separation}{separation checks: \code{"fe"} drops observations in fixed-effect groups
    whose outcome is all zero (iterated with the singleton drop), \code{"relu"} runs the ReLU
    test of Correia, Guimaraes and Zylkin (2019) on the regressors and fixed effects,
    \code{"mu"} flags observations with \eqn{y = 0} whose fitted mean collapses during the
    iterations. \code{"none"} switches all off, \code{"all"} switches all on.}
  \item{keep_singletons}{logical: keep singleton groups. As in Stata, this also switches off
    the fixed-effect separation drop.}
  \item{verbose}{\code{-1} silent, \code{0} notes only, \code{1} or more prints the iteration log.}
}
\details{
Observations are dropped, in this order: missing values in any variable used, zero weights,
singleton groups and groups whose outcome is all zero (iterated over the fixed effects until
nothing changes; both are counted in \code{num_singletons}, as in Stata's ivppmlhdfe with
reghdfe), and observations separated by the ReLU test (counted in \code{num_separated}).
Regressors and instruments that are collinear with the fixed effects or with other columns of
their block are dropped with a note.

Standard errors: the robust variance is
\eqn{(N/(N-1)) B M B} with \eqn{B = (\hat X' W X)^{-1}}; one-way clustering uses
\eqn{G/(G-1)}; multi-way clustering applies the inclusion-exclusion formula with each term's
own \eqn{G/(G-1)} and shifts the matrix by its smallest eigenvalue if that is negative. The
constant is partialled out together with the fixed effects; it is returned in
\code{constant} without a standard error.
}
\value{
An object of class \code{"ivppmlhdfe"}, a list with, among others:
\item{coefficients, vcov, se}{estimates (without the constant), variance matrix and standard errors.}
\item{constant}{the constant (weighted mean of the linear index net of the regressors).}
\item{N}{number of observations used (sum of the weights with frequency weights).}
\item{num_singletons, num_separated}{observations dropped as singletons or separated by a fixed
  effect, and by the ReLU (or mu) separation checks.}
\item{iterations, converged, deviance, ll, ll_0, r2_p}{IRLS iterations, convergence flag,
  deviance, log pseudo-likelihood, null log pseudo-likelihood and pseudo R-squared.}
\item{N_clust}{number of clusters per cluster variable.}
\item{ffirst}{\code{NULL}, or the first-stage statistics: a list with \code{first} (one entry per
  endogenous regressor with F, df1, df2, p, partialR2, SWF, SWdf1, SWdf2, SWp) and the
  system statistics \code{idstat}, \code{iddf}, \code{idp} (Kleibergen-Paap rk LM),
  \code{cdf} (Cragg-Donald Wald F) and \code{widstat} (Kleibergen-Paap rk Wald F).}
\item{esample}{logical vector: rows of \code{data} in the estimation sample.}
\item{eta, mu, d}{linear index, fitted mean and sum of the fixed effects, for the estimation sample.}
}
\references{
Correia, S., Guimaraes, P. and Zylkin, T. (2020). Fast Poisson estimation with high-dimensional
fixed effects. \emph{Stata Journal} 20(1), 95--115.

Correia, S., Guimaraes, P. and Zylkin, T. (2019). Verifying the existence of maximum likelihood
estimates for generalized linear models. arXiv:1903.01633.

Windmeijer, F. and Santos Silva, J. M. C. (1997). Endogeneity in count data models: an
application to demand for health care. \emph{Journal of Applied Econometrics} 12(3), 281--294.
}
\seealso{\code{\link{predict.ivppmlhdfe}}, \code{\link{summary.ivppmlhdfe}}}
\examples{
set.seed(1)
n <- 20; nT <- 8
d <- expand.grid(year = seq_len(nT), id = seq_len(n))
a <- 0.3 * rnorm(n)[d$id] + 0.2 * rnorm(nT)[d$year]
d$z <- rnorm(nrow(d))
d$x2 <- 0.5 * rnorm(nrow(d))
e <- rnorm(nrow(d))
d$x1 <- 0.8 * d$z + e + 0.4 * a
d$y <- exp(0.5 * d$x1 + 0.3 * d$x2 + a + 0.5 * e + 0.8 * rnorm(nrow(d)) - 0.6)

# IV-PPML with two fixed effects and robust standard errors
fit <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d)
summary(fit)

# clustered by id, accel off (reference path)
fit2 <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d, vcov = ~id, accel = FALSE)
coef(fit2)

# first-stage statistics (F, partial R2, Sanderson-Windmeijer F, Kleibergen-Paap, Cragg-Donald)
fit_ff <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d, ffirst = TRUE)
fit_ff$ffirst$first_matrix
c(fit_ff$ffirst$idstat, fit_ff$ffirst$cdf, fit_ff$ffirst$widstat)

# plain PPML (no instruments)
fit3 <- ivppmlhdfe(y ~ x1 + x2 | id + year, data = d)
head(predict(fit3, type = "mu"))
}
\keyword{models}
\keyword{regression}
