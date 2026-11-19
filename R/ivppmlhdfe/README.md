# ivppmlhdfe (R)

IV-PPML with high-dimensional fixed effects, in pure R (no compiled code; base R and `stats` only).
The R port of the Stata command `ivppmlhdfe`, version 0.10.0.

The estimator solves the additive moment E[q(y − μ)] = 0 with q = (x′, z′)′ and
μ = exp(x′β + ΣFE + offset) by iteratively reweighted two-stage least squares, with the fixed
effects concentrated out at each iteration. Without instruments it is plain PPML (as `ppmlhdfe`).

## Install

```r
# install.packages("remotes")
remotes::install_github("ekwonomist/ivppmlhdfe", subdir = "R/ivppmlhdfe")
```

## Example

```r
library(ivppmlhdfe)
set.seed(1)
d <- expand.grid(year = 1:8, id = 1:20)
a <- 0.3 * rnorm(20)[d$id] + 0.2 * rnorm(8)[d$year]
d$z  <- rnorm(nrow(d)); d$x2 <- 0.5 * rnorm(nrow(d)); e <- rnorm(nrow(d))
d$x1 <- 0.8 * d$z + e + 0.4 * a
d$y  <- exp(0.5 * d$x1 + 0.3 * d$x2 + a + 0.5 * e + 0.8 * rnorm(nrow(d)) - 0.6)

# endogenous x1 instrumented by z, id and year fixed effects, clustered by id
fit <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d, vcov = ~id)
summary(fit)

# gravity: exporter-year, importer-year and pair fixed effects, two-way clustering
# ivppmlhdfe(trade ~ rta | exp^year + imp^year + exp^imp | tariff ~ iv, data = g, vcov = ~exp + imp)
```

## Formula

`y ~ exog | fe | endog ~ instruments` (fixest style). Fixed-effect interactions are written
`exp^year`. `y ~ exog | fe` is PPML; `y ~ exog | endog ~ instruments` is IV without fixed effects.

## Main options

| argument | default | meaning |
|---|---|---|
| `vcov` | `"robust"` | or `~pair` (cluster), `~exp + imp` (multi-way cluster) |
| `weights`, `weight_type` | `NULL`, `"pweight"` | weight variable; `"fweight"` for frequency weights |
| `offset`, `exposure` | `NULL` | offset variable, or exposure (offset = log exposure) |
| `accel` | `TRUE` | fast path: warm starts, staged inner tolerance, SE refresh at the converged weights |
| `ffirst` | `FALSE` | first-stage F, partial R², Sanderson-Windmeijer F, Kleibergen-Paap and Cragg-Donald statistics |
| `separation` | `c("fe", "relu")` | separation checks (`"mu"` adds the in-loop check) |
| `tol`, `itol`, `maxit` | `1e-8`, `NULL`, `10000` | IRLS tolerance, inner tolerance, iteration cap |

## Engine

The fixed effect with the most levels is eliminated exactly; the remaining system
S = D₂′W M₁ D₂ is solved by Jacobi-preconditioned conjugate gradients for all columns at once
(the approach of fereg, Schmieder 2026). Group sums are vectorised (`.colSums` on reshaped
balanced groups, `rowsum` otherwise); there are no loops over observations.

## Conventions

Sample drops, separation counts, collinearity checks and standard-error formulas follow the
Stata command (see `?ivppmlhdfe`). On the shared cross-language test set the R results agree
with Stata's ivppmlhdfe v0.9.4 to better than 1e-7 (coefficients and standard errors), with
identical observation, singleton and separation counts.

License: MIT. Authors: Ohyun Kwon, Mario Larch, Jangsu Yoon, Yoto V. Yotov.
