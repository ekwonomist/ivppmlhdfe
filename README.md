> **GitHub version.** For the local browsing version with full annotated code examples, open [`readme.html`](readme.html).

# ivppmlhdfe

Instrumental-variables Poisson pseudo-maximum likelihood estimation with high-dimensional fixed effects.

**Alpha release** (v0.9.4) — shared with coauthors for testing. Please report issues to Ohyun Kwon.

Ohyun Kwon, Mario Larch, Jangsu Yoon, Yoto V. Yotov

## Overview

`ivppmlhdfe` estimates IV-PPML models with multiple sets of high-dimensional fixed effects. The estimator targets the additive moment condition E[Z(y - μ)] = 0 ([Mullahy, 1997](https://doi.org/10.2307/2951380), Eq. 6), solved via iteratively reweighted 2SLS with fixed effects concentrated out at each iteration using [`reghdfe`](https://github.com/sergiocorreia/reghdfe).

The command is designed to feel natural to users of [`ppmlhdfe`](https://github.com/sergiocorreia/ppmlhdfe). All standard post-estimation commands are supported.

| Class | Dimension | FE Structure | Example |
|-------|-----------|--------------|---------|
| A | N × T | Individual + Time | `absorb(id year)` |
| B | N × N × T | Exporter-year + Importer-year | `absorb(exp#yr imp#yr)` |
| C | N × N × T | Pair + Exporter-year + Importer-year | `absorb(exp#imp exp#yr imp#yr)` |

**Note:** Both **just-identified** and **overidentified** models are supported. Multiple endogenous regressors are allowed.

## Installation

1. [Download the latest release](https://github.com/ekwonomist/ivppmlhdfe/archive/refs/heads/main.zip) and extract to a local folder.

2. Install dependencies:
```stata
ssc install reghdfe, replace
ssc install ppmlhdfe, replace
ssc install ftools, replace
```

3. Add the folder to your Stata ado path:
```stata
adopath ++ "path/to/ivppmlhdfe"
```

## Syntax

```stata
ivppmlhdfe depvar [exogvars] (endogvar = instrument) [if] [in] [pw/fw], absorb(absvars) [options]
```

### Options

| Option | Description |
|--------|-------------|
| `absorb(absvars)` | Fixed effects to absorb (required). Factor syntax supported. |
| `noabsorb` | Estimate without fixed effects (mutually exclusive with `absorb()`). |
| `vce(robust\|cluster c1 [c2 [c3]])` | Variance estimator. Default: robust. Multi-way clustering supported. |
| `exposure(varname)` | Exposure variable (log used as offset). |
| `offset(varname)` | Offset variable with coefficient fixed at 1. |
| `d(newvarname)` / `d2` | Save sum of absorbed FE. `d2` auto-names the variable. Required for `predict mu`. |
| `separation(default\|all\|none\|fe simplex relu mu)` | Separation detection. Default: `fe simplex relu`. `all` adds in-loop `mu` detection. |
| `standardize` | Standardize X and Z columns to unit variance for numerical stability (results identical after back-transform). |
| `guess(simple\|mean)` | Initial mu guess. `simple` (default) = 0.5*(y + mean(y)); `mean` = mean(y). |
| `keepsingletons` | Keep singleton observations (default drops them). |
| `tagsep(newvarname)` | Tag separated obs and exit (no estimation). |
| `zvarname(newvarname)` | Save ReLU separation certificate variable. |
| `eform` / `irr` | Display exponentiated coefficients (incidence-rate ratios). |
| `tolerance(#)` | Outer convergence tolerance. Default: 1e-8. |
| `itolerance(#)` | Inner FE-solver target tolerance. Default: max(1e-12, 0.1*tol). |
| `maxiterations(#)` | Max IRLS iterations. Default: 10000. |
| `nolog` | Suppress IRLS iteration log. |
| `verbose(#)` | Verbosity: -1 (silent), 0 (default), 1 (iterations), 2 (detail). |

Unrecognised options (e.g. `dof()`, `pool()`, `accel()`) are forwarded to `reghdfe` via the `absorb()` string.

### Predict

After estimation, `predict` supports:

| Statistic | Formula | Requires `d()`? |
|-----------|---------|-----------------|
| `mu` (default) | exp(xb + d) | Yes |
| `xb` | X*β + _cons + offset | No |
| `xbd` / `eta` | xb + d | Yes |
| `d` | Sum of fixed effects | Yes |
| `scores` / `residuals` | y - μ | Yes |
| `pearson` | (y - μ) / √μ | Yes |
| `working` | (y - μ) / μ | Yes |
| `deviance` | Deviance residual | Yes |

### Post-Estimation Commands

All verified: `test`, `testparm`, `lincom`, `nlcom`, `predict`, `margins`, `estat ic`, `estat summarize`, `estat vce`, `regsave`, `esttab`/`estout`, `outreg2`, `coefplot`.

## Feature Comparison with ppmlhdfe

| # | Feature | ppmlhdfe | ivppmlhdfe |
|---|---------|:--------:|:----------:|
| 1 | Multi-way clustering | ✓ | ✓ |
| 2 | Separation detection (fe) | ✓ | ✓ |
| 3 | Singleton drop | ✓ | ✓ |
| 4 | predict (mu, eta, scores, residuals, ...) | ✓ | ✓ |
| 5 | margins | ✓ | ✓ |
| 6 | estat ic / summarize / vce | ✓ | ✓ |
| 7 | exposure() / offset() | ✓ | ✓ |
| 8 | d() option | ✓ | ✓ |
| 9 | eform / irr | ✓ | ✓ |
| 10 | Factor variable syntax in absorb() | ✓ | ✓ |

## Examples

```stata
* Two-way gravity model
ivppmlhdfe trade (policy = instrument), absorb(exp_year imp_year) vce(cluster pair_id)

* Three-way gravity model with pair FE
ivppmlhdfe trade (policy = instrument), absorb(exp_year imp_year pair_id) vce(cluster pair_id)

* With exposure and prediction
ivppmlhdfe trade (policy = instrument), absorb(exp_year imp_year) exposure(pop) d(fe_sum)
predict muhat, mu
predict resid, residuals

* Multi-way clustering
ivppmlhdfe trade (policy = instrument), absorb(exp_year imp_year) vce(cluster exp imp)

* Incidence-rate ratios
ivppmlhdfe trade (policy = instrument), absorb(exp_year imp_year) vce(cluster pair_id) irr

* Margins
ivppmlhdfe y x2 (x1 = z), absorb(id year) d(fe_sum)
margins, dydx(x2)
```

## Julia Backend (`ivppmlhdfejl`) — Optional

A Julia-powered backend is available for faster estimation. It uses [`FixedEffects.jl`](https://github.com/FixedEffects/FixedEffects.jl) for FE absorption with solver reuse across IRLS iterations.

**Pure Julia is ~7x faster than Stata/Mata** for repeated estimation (Monte Carlo, bootstrap).

### Requirements

- [Julia](https://julialang.org/downloads/) 1.9 or later
- Stata 15.0 or later
- [`jl`](https://github.com/droodman/julia.ado) (Stata-Julia bridge by David Roodman)

### Installation

1. Install Julia (1.9+) and the `jl` Stata-Julia bridge:
```stata
ssc install julia, replace
```

2. Run the automated installer (copies bridge files + Julia module to your ado directory):
```stata
do "C:/your_path/ivppmlhdfe/julia/install_ivppmlhdfejl.do"
```

3. First run will install Julia dependencies and precompile (~30s one-time cost).

### Usage

```stata
* IV Poisson PML (same syntax as ivppmlhdfe)
ivppmlhdfejl trade (policy = instrument), absorb(exp_year imp_year) vce(cluster pair_id)

* Exogenous Poisson PML (no endogenous variable — runs ppmlhdfejl)
ivppmlhdfejl trade x1 x2, absorb(exp_year imp_year) vce(cluster pair_id)
```

### Three Ways to Run

| Method | Command | Speed | Best for |
|--------|---------|-------|----------|
| Stata/Mata | `ivppmlhdfe` | 1x | General use |
| Julia via Stata | `ivppmlhdfejl` | Varies | Large single datasets |
| Pure Julia | `IVPPMLFixedEffectModels.jl` | ~7x | Monte Carlo, bootstrap |

## SPJ Bias Correction + Bootstrap SE

IV-PPML with fixed effects suffers from incidental parameter bias that does not arise in standard PPML. We recommend split-panel jackknife (SPJ) bias correction paired with bootstrap standard errors.

Ready-to-run do-files with example data are included in the `data/` folder:

| Class | SPJ Method | Formula | Do-file |
|-------|-----------|---------|---------|
| A | FVW 2016 (time + individual split) | beta_SPJ = 3*b_full - b_Thalf - b_Nhalf | [`MC_SPJ_BTS_ClassA.do`](data/MC_SPJ_BTS_ClassA.do) |
| B | WZ 2021 (4-subpanel country split) | beta_SPJ = 2*b_full - mean(b_aa, b_ab, b_ba, b_bb) | [`MC_SPJ_BTS_ClassB.do`](data/MC_SPJ_BTS_ClassB.do) |
| C | 8-panel (country x time split) | beta_8p = 4*b_full - 2*b_country - 2*b_time + b_8cell | [`MC_SPJ_BTS_ClassC.do`](data/MC_SPJ_BTS_ClassC.do) |

Each script loads example data (`ivppmlhdfe_ClassX.dta`), computes the SPJ point estimate, and runs B=200 bootstrap replications to produce SE and 95% CI. See [`readme.html`](readme.html) for full annotated code with step-by-step explanations.

To generate fresh example data, run [`data/DGP.do`](data/DGP.do).

**Computational cost:** SPJ requires 3-9 estimations per dataset, times B bootstrap replications. With B=200, expect 600-1,800 `ivppmlhdfe` calls per dataset.

## Example Data

The `data/` folder contains:

| File | Description |
|------|-------------|
| `DGP.do` | Data generating process (creates all three datasets below) |
| `ivppmlhdfe_ClassA.dta` | N=100, T=10 panel with individual + time FE |
| `ivppmlhdfe_ClassB.dta` | Nc=30, T=20 gravity panel with two-way FE |
| `ivppmlhdfe_ClassC.dta` | Nc=30, T=20 gravity panel with three-way FE |
| `MC_SPJ_BTS_ClassA.do` | SPJ + bootstrap template for Class A |
| `MC_SPJ_BTS_ClassB.do` | SPJ + bootstrap template for Class B |
| `MC_SPJ_BTS_ClassC.do` | SPJ + bootstrap template for Class C |

## Known Limitations

- **Separation detection**: All three separation detection methods (fe, simplex, relu) are supported, matching ppmlhdfe. The default applies all three in sequence.
- **Convergence failures**: rc=430 (non-convergence) and rc=9003 (runaway divergence) are now loud failures with clear error messages, matching Stata's native ivppmlhdfe.
- **Collinearity detection**: v0.9.4 automatically detects and drops collinear regressors/instruments (including factor-variable dummies fully absorbed by FE).
- **Constant term**: `_cons` is reported without a standard error (by design; use bootstrap for `_cons` inference).
- **No first-stage diagnostics**: No Cragg-Donald or Kleibergen-Paap statistics.
- **Incidental parameter bias**: IV-PPML has O(1/T) + O(1/N) point-estimate bias in short panels. Use the SPJ + bootstrap templates above.

## References

- Mullahy, J. (1997). "Instrumental-variable estimation of count data models." *Review of Economics and Statistics*, 79(4), 586-593.
- Correia, S., P. Guimaraes, and T. Zylkin (2020). "Fast Poisson estimation with high-dimensional fixed effects." *Stata Journal*, 20(1), 95-115.
- Weidner, M. and T. Zylkin (2021). "Bias and consistency in three-way gravity models." *Journal of International Economics*, 132, 103513.
- Fernandez-Val, I. and M. Weidner (2016). "Individual and time effects in nonlinear panel models with large N, T." *Journal of Econometrics*, 192(1), 291-312.
- Cameron, A. C., J. B. Gelbach, and D. L. Miller (2011). "Robust inference with multiway clustering." *Journal of Business and Economic Statistics*, 29(2), 238-249.

## Authors

Ohyun Kwon, Mario Larch, Jangsu Yoon, Yoto V. Yotov

Bug reports: [GitHub Issues](https://github.com/ekwonomist/ivppmlhdfe/issues)
