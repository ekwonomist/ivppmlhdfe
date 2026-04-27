# ivppmlhdfe

Instrumental-variable Poisson pseudo-maximum likelihood estimation with high-dimensional fixed effects.

v0.9.4 — please report issues to Ohyun Kwon at theekwonomist@gmail.com.

Ohyun Kwon, Mario Larch, Jangsu Yoon, Yoto V. Yotov

## Overview

`ivppmlhdfe` estimates IV-PPML models with multiple sets of high-dimensional fixed effects. The estimator targets the additive moment condition E[q(y - μ)] = 0 with q = (x', z')' (exogenous regressors stacked with excluded instruments), together with the per-group fixed-effect score Σ_{g∈r}(y_g - μ_g) = 0, following [Windmeijer and Santos Silva (1997)](https://doi.org/10.1002/(SICI)1099-1255(199705)12:3%3C281::AID-JAE436%3E3.0.CO;2-1). It is solved via iteratively reweighted 2SLS [(Correia, Guimaraes, and Zylkin, 2020)](https://journals.sagepub.com/doi/10.1177/1536867X20909691) with fixed effects concentrated out at each iteration using [`reghdfe`](https://github.com/sergiocorreia/reghdfe).

For incidental-parameter bias under IV-PPML and the split-panel jackknife (SPJ) + bootstrap remedy, see the companion paper [Kwon, Larch, Yoon, and Yotov (2026)](https://ideas.repec.org/p/drx/wpaper/202611.html).

The command is designed to feel natural to users of [`ppmlhdfe`](https://github.com/sergiocorreia/ppmlhdfe). All standard post-estimation commands are supported. Both **just-identified** and **overidentified** models are supported, and multiple endogenous regressors are allowed.

| Class | Dimension | FE Structure | Example |
|-------|-----------|--------------|---------|
| A | N × T | Individual + Time | `absorb(id year)` |
| B | N × N × T | Exporter-year + Importer-year | `absorb(exp#year imp#year)` |
| C | N × N × T | Pair + Exporter-year + Importer-year | `absorb(exp#imp exp#year imp#year)` |

## Installation

Install dependencies:
```stata
ssc install reghdfe, replace
ssc install ppmlhdfe, replace
ssc install ftools, replace
```

Install `ivppmlhdfe`:
```stata
net install ivppmlhdfe, from("https://raw.githubusercontent.com/ekwonomist/ivppmlhdfe/main/") replace
```

> **Example data.** The `data/` folder (example datasets and SPJ+bootstrap templates) is not bundled with `net install`. To use the examples, clone or download the repository manually from [GitHub](https://github.com/ekwonomist/ivppmlhdfe).

## Quick Example

```stata
* Individual + time fixed effects (Class A)
use "data/ivppmlhdfe_ClassA.dta", clear
ivppmlhdfe y x2 (x1 = z), absorb(id year) vce(robust)

* Two-way gravity (Class B)
use "data/ivppmlhdfe_ClassB.dta", clear
ivppmlhdfe y x2 (x1 = z), absorb(exp#year imp#year) vce(cluster pair)

* Three-way gravity with pair FE (Class C)
use "data/ivppmlhdfe_ClassC.dta", clear
ivppmlhdfe y x2 (x1 = z), absorb(exp#imp exp#year imp#year) vce(cluster pair)
```

For full option documentation, run `help ivppmlhdfe` after installation.

## Syntax

```stata
ivppmlhdfe depvar [exogvars] (endogvars = instruments) [if] [in] [pw/fw], absorb(absvars) [options]
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

## Further Examples

```stata
* With exposure and prediction
ivppmlhdfe y x2 (x1 = z), absorb(exp#year imp#year) exposure(pop) d(fe_sum)
predict muhat, mu
predict resid, residuals

* Multi-way clustering
ivppmlhdfe y x2 (x1 = z), absorb(exp#year imp#year) vce(cluster exp imp)

* Incidence-rate ratios
ivppmlhdfe y x2 (x1 = z), absorb(exp#year imp#year) vce(cluster pair) irr

* Margins
ivppmlhdfe y x2 (x1 = z), absorb(id year) d(fe_sum)
margins, dydx(x2)
```

## Julia Backend (`ivppmlhdfejl`) — Optional

A Julia-powered backend is available for faster estimation. It uses [`FixedEffects.jl`](https://github.com/FixedEffects/FixedEffects.jl) for FE absorption with solver reuse across IRLS iterations.

Pure Julia is substantially faster than Stata/Mata for repeated estimation (Monte Carlo, bootstrap), since it avoids the per-call overhead of the Stata-Julia bridge.

### Requirements

- [Julia](https://julialang.org/downloads/) 1.9 or later
- Stata 15.0 or later
- [`jl`](https://github.com/droodman/julia.ado) (Stata-Julia bridge by David Roodman)

### Installation

1. Install Julia (1.9+) and the `jl` Stata-Julia bridge:
```stata
ssc install julia, replace
```

2. From your cloned/downloaded ivppmlhdfe repo, run the installer:
```stata
do "julia/install_ivppmlhdfejl.do"
```

3. First run will install Julia dependencies and precompile (~30s one-time cost).

### Usage

```stata
* IV Poisson PML (same syntax as ivppmlhdfe)
ivppmlhdfejl y x2 (x1 = z), absorb(exp#year imp#year) vce(cluster pair)

* Exogenous Poisson PML (no endogenous variable — runs ppmlhdfejl)
ivppmlhdfejl y x2 x1, absorb(exp#year imp#year) vce(cluster pair)
```

### Three Ways to Run

| Method | Command | Speed | Best for |
|--------|---------|-------|----------|
| Stata/Mata | `ivppmlhdfe` | 1x | General use |
| Julia via Stata | `ivppmlhdfejl` | Varies | Large single datasets |
| Pure Julia | `IVPPMLFixedEffectModels.jl` | Fastest | Monte Carlo, bootstrap |

## SPJ Bias Correction + Bootstrap SE

IV-PPML with fixed effects suffers from incidental parameter bias that does not arise in standard PPML. The companion paper [Kwon, Larch, Yoon, and Yotov (2026)](https://ideas.repec.org/p/drx/wpaper/202611.html) derives the bias orders by FE structure and develops a split-panel jackknife (SPJ) bias correction paired with bootstrap standard errors. See the paper for the full derivations, the recommended bootstrap aggregator (median CI-implied SE), and the Monte Carlo evidence.

Ready-to-run do-files with example data are included in the `data/` folder (clone the repo from GitHub to access them):

| Class | SPJ Method | Formula | Do-file |
|-------|-----------|---------|---------|
| A | FVW 2016 (time + individual split) | beta_SPJ = 3*b_full - b_Thalf - b_Nhalf | [`MC_SPJ_BTS_ClassA.do`](data/MC_SPJ_BTS_ClassA.do) |
| B | WZ 2021 (4-subpanel country split) | beta_SPJ = 2*b_full - mean(b_aa, b_ab, b_ba, b_bb) | [`MC_SPJ_BTS_ClassB.do`](data/MC_SPJ_BTS_ClassB.do) |
| C | 8-panel (country × time split) | beta_8p = 4*b_full - 2*b_country - 2*b_time + b_8cell | [`MC_SPJ_BTS_ClassC.do`](data/MC_SPJ_BTS_ClassC.do) |

Each script loads example data (`ivppmlhdfe_ClassX.dta`), computes the SPJ point estimate, and runs B=1000 bootstrap replications (matching the paper's MC design) to produce SE and 95% CI.

To generate fresh example data, run [`data/DGP.do`](data/DGP.do).

**Computational cost:** SPJ requires 3–8 estimations per dataset, times B bootstrap replications. With B=1000, expect 3,000–8,000 `ivppmlhdfe` calls per dataset.

## Example Data

The `data/` folder (available when you clone the repository) contains:

| File | Description |
|------|-------------|
| `DGP.do` | Data generating process (creates all three datasets below) |
| `ivppmlhdfe_ClassA.dta` | N=100, T=10 panel with individual + time FE |
| `ivppmlhdfe_ClassB.dta` | Nc=30, T=20 gravity panel with two-way FE |
| `ivppmlhdfe_ClassC.dta` | Nc=30, T=20 gravity panel with three-way FE |
| `MC_SPJ_BTS_ClassA.do` | SPJ + bootstrap template for Class A |
| `MC_SPJ_BTS_ClassB.do` | SPJ + bootstrap template for Class B |
| `MC_SPJ_BTS_ClassC.do` | SPJ + bootstrap template for Class C |

## References

- Cameron, A. C., J. B. Gelbach, and D. L. Miller (2011). "Robust inference with multiway clustering." *Journal of Business and Economic Statistics*, 29(2), 238–249.
- Correia, S., P. Guimaraes, and T. Zylkin (2020). "Fast Poisson estimation with high-dimensional fixed effects." *Stata Journal*, 20(1), 95–115.
- Fernandez-Val, I. and M. Weidner (2016). "Individual and time effects in nonlinear panel models with large N, T." *Journal of Econometrics*, 192(1), 291–312.
- Kwon, O., M. Larch, J. Yoon, and Y. V. Yotov (2026). "Instrumental-Variable Poisson PML with High-Dimensional Fixed Effects." *Drexel University School of Economics Working Paper* 2026-11. [[ideas.repec.org]](https://ideas.repec.org/p/drx/wpaper/202611.html)
- Mullahy, J. (1997). "Instrumental-variable estimation of count data models." *Review of Economics and Statistics*, 79(4), 586–593. *(Related but distinct: a transformation estimator for multiplicative unobserved heterogeneity, not the additive moment used here.)*
- Weidner, M. and T. Zylkin (2021). "Bias and consistency in three-way gravity models." *Journal of International Economics*, 132, 103513.
- Windmeijer, F. A. G., and J. M. C. Santos Silva (1997). "Endogeneity in count data models: An application to demand for health care." *Journal of Applied Econometrics*, 12(3), 281–294.

Report issues to Ohyun Kwon (theekwonomist@gmail.com) or open a [GitHub issue](https://github.com/ekwonomist/ivppmlhdfe/issues).
