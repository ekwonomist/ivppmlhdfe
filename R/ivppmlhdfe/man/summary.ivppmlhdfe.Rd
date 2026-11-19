\name{summary.ivppmlhdfe}
\alias{summary.ivppmlhdfe}
\alias{print.ivppmlhdfe}
\alias{print.summary.ivppmlhdfe}
\alias{coef.ivppmlhdfe}
\alias{vcov.ivppmlhdfe}
\alias{nobs.ivppmlhdfe}
\alias{deviance.ivppmlhdfe}
\alias{logLik.ivppmlhdfe}
\alias{fitted.ivppmlhdfe}
\title{Methods for ivppmlhdfe fits}
\description{Coefficient table, coefficients, variance matrix and other accessors.}
\usage{
\method{summary}{ivppmlhdfe}(object, level = 0.95, ...)
\method{print}{summary.ivppmlhdfe}(x, digits = max(3L, getOption("digits") - 3L), ...)
\method{print}{ivppmlhdfe}(x, digits = max(3L, getOption("digits") - 3L), ...)
\method{coef}{ivppmlhdfe}(object, ...)
\method{vcov}{ivppmlhdfe}(object, ...)
\method{nobs}{ivppmlhdfe}(object, ...)
\method{deviance}{ivppmlhdfe}(object, ...)
\method{logLik}{ivppmlhdfe}(object, ...)
\method{fitted}{ivppmlhdfe}(object, ...)
}
\arguments{
  \item{object, x}{an \code{"ivppmlhdfe"} fit (or its summary).}
  \item{level}{confidence level of the intervals in the coefficient table.}
  \item{digits}{significant digits to print.}
  \item{...}{unused.}
}
\details{
Inference uses the normal distribution (z statistics), as the Stata command. \code{coef} and
\code{vcov} exclude the partialled-out constant; it is stored in \code{object$constant}.
\code{fitted} returns the fitted means (\code{NA} outside the estimation sample).
}
\value{
\code{summary} returns the fit with a \code{coeftable} element (estimate, standard error, z,
p-value and confidence interval); the accessors return the corresponding values.
}
\examples{
set.seed(3)
d <- expand.grid(year = 1:6, id = 1:15)
d$x <- rnorm(nrow(d))
d$y <- rpois(nrow(d), exp(0.5 * d$x + rnorm(15)[d$id]))
fit <- ivppmlhdfe(y ~ x | id + year, data = d, vcov = ~id)
summary(fit)
coef(fit); vcov(fit); nobs(fit); logLik(fit)
}
