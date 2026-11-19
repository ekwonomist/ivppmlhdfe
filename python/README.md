# ivppmlhdfe (Python)

IV-PPML: instrumental-variables Poisson pseudo-maximum likelihood with high-dimensional fixed
effects. A native Python port of the Stata command `ivppmlhdfe` (same estimator, same conventions
for starting values, convergence, singletons, separation, collinearity and standard errors), with
its own fixed-effect engine. Dependencies: numpy, scipy, pandas.

## Install

```
pip install "git+https://github.com/ekwonomist/ivppmlhdfe#subdirectory=python"
```

## Example

```python
import pandas as pd
import ivppmlhdfe

df = pd.read_csv("trade.csv")
r = ivppmlhdfe.fit(df, "y ~ x2 | exp^year + imp^year + exp^imp | x1 ~ z",
                   vcov={"CL": ["pair"]})
print(r.summary())
r.coef, r.se, r.vcov          # pandas Series / DataFrame
mu = r.predict("mu")          # in-sample fitted values (NaN outside the estimation sample)
```

Formula grammar (fixest style): `depvar ~ exogenous | fixed effects | endogenous ~ instruments`.
Fixed-effect interactions are written `exp^year`. The FE part and the IV part are optional: no
IV part is plain PPML (the same estimates as Stata's `ppmlhdfe`).

Main options of `ivppmlhdfe.fit(data, formula, ...)`:

| option | default | meaning |
|---|---|---|
| `vcov` | `"robust"` | `"robust"` or `{"CL": ["pair"]}`; two- and three-way clustering with several names |
| `weights`, `weight_type` | `None`, `"pweight"` | weight column; `"pweight"`, `"aweight"` or `"fweight"` |
| `offset`, `exposure` | `None` | offset column, or exposure column (offset = log(exposure)) |
| `accel` | `True` | fast path (warm start, fereg tolerance schedule, VCE refresh at the converged weights); `False` gives the v0.9.4 path |
| `ffirst` | `False` | first-stage statistics (F, partial R², Sanderson-Windmeijer F, Kleibergen-Paap, Cragg-Donald) |
| `tol`, `itol`, `maxit` | `1e-8`, `max(1e-12, 0.1 tol)`, `10000` | IRLS tolerance on the deviance, inner FE tolerance, iteration cap |
| `separation` | `("fe", "relu")` | separation checks: `"fe"`, `"relu"`, `"mu"` (or `"all"`, `"none"`) |
| `keep_singletons` | `False` | keep singleton groups |

Result fields: `coef`, `se`, `vcov`, `intercept` (the absorbed constant, reported without a
standard error), `N`, `N_full`, `iterations`, `converged`, `deviance`, `ll`, `ll_0`,
`num_singletons` (singletons and observations in all-zero FE groups, as in Stata), `num_separated`
(ReLU and mu separation), `N_clust`, `ffirst`, `stats` (timings, inner iterations).

With `ffirst=True`, `r.ffirst` holds, for each endogenous regressor, the first-stage
`F`, `df1`, `df2`, `p`, `partialR2` and the Sanderson-Windmeijer `SWF`, `SWdf1`, `SWdf2`, `SWp`,
and for the system the Kleibergen-Paap rk LM statistic `idstat` (`iddf`, `idp`), the Cragg-Donald
Wald F `cdf` and the Kleibergen-Paap rk Wald F `widstat`, with the conventions of
`ivreghdfe ..., ffirst` at the final IRLS weights.

## Fixed-effect engine

The fixed effect with the most levels is eliminated exactly (weighted within-group demeaning);
the remaining fixed effects are solved by Jacobi-preconditioned conjugate gradients on the Schur
complement, all columns at once, with fereg's certificate check (Schmieder 2026, NBER WP 36017).
With one fixed effect this is the exact within transform. All loops are vectorised in numpy and
scipy.sparse.

## Authors and license

Ohyun Kwon, Mario Larch, Jangsu Yoon, Yoto V. Yotov. MIT license.
