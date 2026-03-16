# ivppmlhdfe

Instrumental-variables Poisson pseudo-maximum likelihood estimation with high-dimensional fixed effects.

Ohyun Kwon, Mario Larch, Jangsu Yoon, Yoto V. Yotov

## Overview

`ivppmlhdfe` estimates IV-PPML models with multiple sets of high-dimensional fixed effects using iteratively reweighted GMM ([Mullahy, 1997](https://doi.org/10.2307/2951380)). Fixed effects are concentrated out via the PPML first-order condition at each iteration using [`reghdfe`](https://github.com/sergiocorreia/reghdfe).

The command handles panel models applied in applied econometrics:
| Class | Dimension | FE Structure |
|-------|-----------|--------------|
| RE | N x T | Time |
| A | N x T | Individual + Time |
| B | N x N x T | Exporter-year + Importer-year |
| C | N x N x T | Exporter-year + Importer-year + Pair |

## Installation

> ![Note](https://img.shields.io/badge/Note-Private_git_limitation-orange) `net install` is not available.

1. [Download the latest release](https://github.com/ekwonomist/ivppmlhdfe/archive/refs/heads/main.zip) and copy `ivppmlhdfe.ado` and `ivppmlhdfe.sthlp` to your Stata ado directory (e.g., type `sysdir` in Stata to find it).

2. Install dependencies:
```stata
ssc install reghdfe, replace
ssc install ppmlhdfe, replace
ssc install ftools, replace
```

## Syntax

```stata
ivppmlhdfe depvar [exogvars] (endogvars = instruments) [if] [in] [pw], absorb(absvars) [options]
```

### Options

| Option | Description |
|--------|-------------|
| `absorb(absvars)` | Fixed effects to absorb (required) |
| `vce(robust\|cluster clustvar)` | Variance estimator (default: robust) |
| `tolerance(#)` | Convergence tolerance (default: 1e-8) |
| `maxiterations(#)` | Maximum IRLS iterations (default: 1000) |
| `biascorrection(class id1 id2)` | Analytical bias correction (class `a` or `b`) |
| `eform` / `irr` | Display exponentiated coefficients |
| `verbose(#)` | Verbosity level 0-3 |

## Examples

> ![Note](https://img.shields.io/badge/Note-Private_git_limitation-orange) `net install` is not available. Example datasets are included in the `data/` folder.

```stata
* Class RE: Time FE only
use "data/ivppmlhdfe_ClassRE.dta", clear
ivppmlhdfe y (x = z), absorb(year) vce(robust)

* Class A: Individual + Time FE
use "data/ivppmlhdfe_ClassA.dta", clear
ivppmlhdfe y (x = z), absorb(id year) vce(cluster id)

* Class B: Exporter-year + Importer-year FE
use "data/ivppmlhdfe_ClassB.dta", clear
ivppmlhdfe trade (policy = instrument), absorb(exp_year imp_year) vce(cluster pair_id)

* Class C: Exporter-year + Importer-year + Pair FE
use "data/ivppmlhdfe_ClassC.dta", clear
ivppmlhdfe trade (policy = instrument), absorb(exp_year imp_year pair_id) vce(cluster pair_id)
```

## Julia Backend (`ivppmlhdfejl`) — optional

A Julia-powered backend is available for faster estimation on large datasets. It uses [`FixedEffects.jl`](https://github.com/FixedEffects/FixedEffects.jl) for FE absorption with solver reuse across IRLS iterations. This section can be skipped if you only need the Stata backend.

### Requirements

- [Julia](https://julialang.org/downloads/) 1.9 or later
- Stata 15.0 or later
- [`jl`](https://github.com/droodman/julia.ado) (Stata-Julia bridge by David Roodman)

### Installation

> ![Note](https://img.shields.io/badge/Note-Private_git_limitation-orange) Manual file copy is required instead of the one-liner install script.

1. Install Julia and the `jl` Stata-Julia bridge:
```stata
ssc install julia, replace
```

2. Copy `julia/ivppmlhdfejl.ado`, `julia/ivppmlhdfejl_load.ado`, `julia/ivppmlhdfejl_project.toml`, and the `julia/IVPPMLFixedEffectModels/` folder to your Stata ado directory.

3. First run will install Julia dependencies and precompile (one-time cost).

### Usage

```stata
* Same syntax as ivppmlhdfe, just replace the command name
ivppmlhdfejl trade (policy = instrument), absorb(exp_year imp_year) vce(cluster pair_id)
```

### Current Limitations

- No bias correction (`biascorrection()` option not yet implemented)
- No variance correction (Kauermann-Carroll adjustment not available)
- Requires `jl` (Stata-Julia bridge) to be installed and configured
- First invocation in a Stata session has a startup cost (~10-20 seconds for Julia initialization)
- GPU acceleration (`gpu` option) is experimental

### Performance

Benchmarked at approximately 1.8x faster than the Stata/Mata backend for three-way gravity models (N=200 countries, T=50 periods).

## ![Development Status](https://img.shields.io/badge/Development_Status-In_Progress-black)

This package is under active development. The core estimator is functional, but several components are in progress:

- [x] Core IV-PPML estimator (IRLS-GMM with `reghdfe`)
- [x] Robust and cluster-robust variance estimation
- [x] Multi-way fixed effects (Classes A, B, C)
- [x] Kauermann-Carroll variance correction
- [x] Julia backend (`ivppmlhdfejl`)
- [ ] Point-estimate bias correction (split-panel jackknife)
- [ ] Variance correction for jackknife estimates
- [ ] Class C bias correction
- [ ] Stata help file examples with real data
- [ ] SSC submission

### ![Known Issues](https://img.shields.io/badge/Known_Issues-black)

- **Incidental-parameter bias**: IV-PPML has O(1/T) + O(1/N) point-estimate bias in Class A models (individual + time FE) and O(1/N) in Class B (directional FE). Standard PPML does not have this bias due to Bartlett identity cancellation, but IV-PPML breaks this cancellation because instruments are not likelihood scores. Split-panel jackknife correction is in development.
- **Variance undercoverage**: Cluster-robust standard errors can undercover in short panels (coverage ~92% vs nominal 95%). Improved variance correction is in progress.

## References

- Mullahy, J. (1997). "Instrumental-variable estimation of count data models." *Review of Economics and Statistics*, 79(4), 586-593.
- Correia, S., P. Guimaraes, and T. Zylkin (2020). "Fast Poisson estimation with high-dimensional fixed effects." *Stata Journal*, 20(1), 95-115.
- Weidner, M. and T. Zylkin (2021). "Bias and consistency in three-way gravity models." *Journal of International Economics*, 132, 103513.
- Fernandez-Val, I. and M. Weidner (2016). "Individual and time effects in nonlinear panel models with large N, T." *Journal of Econometrics*, 192(1), 291-312.
