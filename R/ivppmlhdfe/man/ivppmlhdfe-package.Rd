\name{ivppmlhdfe-package}
\alias{ivppmlhdfe-package}
\docType{package}
\title{IV-PPML with High-Dimensional Fixed Effects}
\description{
Instrumental-variables Poisson pseudo-maximum likelihood with any number of high-dimensional
fixed effects, in pure R. The estimator solves the additive moment
\eqn{E[q(y - \mu)] = 0}, \eqn{q = (x', z')'}, \eqn{\mu = \exp(x'\beta + \sum FE + offset)}, by
iteratively reweighted two-stage least squares with the fixed effects concentrated out at each
iteration (Windmeijer and Santos Silva 1997; Correia, Guimaraes and Zylkin 2020). It is the R
port of the Stata command \code{ivppmlhdfe} (version 0.10.0 of the shared code base) and follows
its conventions for singletons, separation, collinearity, weights, offsets and standard errors.
The main function is \code{\link{ivppmlhdfe}}.
}
\author{
Ohyun Kwon, Mario Larch, Jangsu Yoon, Yoto V. Yotov. Maintainer: Ohyun Kwon <theekwonomist@gmail.com>.
}
\keyword{package}
