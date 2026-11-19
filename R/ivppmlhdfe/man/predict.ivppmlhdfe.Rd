\name{predict.ivppmlhdfe}
\alias{predict.ivppmlhdfe}
\title{Predictions from an ivppmlhdfe fit}
\description{
Fitted means, linear predictions and the fixed-effect sum, as \code{predict} after the Stata
command (\code{mu}, \code{xb}, \code{eta}, \code{d}).
}
\usage{
\method{predict}{ivppmlhdfe}(object, newdata = NULL, type = c("mu", "xb", "eta", "d"), ...)
}
\arguments{
  \item{object}{an \code{"ivppmlhdfe"} fit.}
  \item{newdata}{\code{NULL} (the estimation data: one value per row of the original data,
    \code{NA} outside the estimation sample) or a data frame. With new data only
    \code{type = "xb"} is available when the model has fixed effects.}
  \item{type}{\code{"mu"} fitted mean \eqn{\exp(\eta)}; \code{"xb"} the regressors' index
    including the constant; \code{"eta"} the full linear index \eqn{xb + d + offset};
    \code{"d"} the sum of the fixed effects (weighted mean zero).}
  \item{...}{unused.}
}
\value{A numeric vector.}
\examples{
set.seed(2)
d <- expand.grid(year = 1:6, id = 1:15)
d$x <- rnorm(nrow(d))
d$y <- rpois(nrow(d), exp(0.5 * d$x + rnorm(15)[d$id]))
fit <- ivppmlhdfe(y ~ x | id + year, data = d)
summary(predict(fit, type = "mu") - d$y)
}
