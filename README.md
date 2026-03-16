# ivppmlhdfe

Instrumental-variables Poisson pseudo-maximum likelihood estimation with high-dimensional fixed effects.

## Overview

`ivppmlhdfe` estimates IV-PPML models with multiple sets of high-dimensional fixed effects using iteratively reweighted GMM ([Mullahy, 1997](https://doi.org/10.2307/2951380)). Fixed effects are concentrated out via the PPML first-order condition at each iteration using [`reghdfe`](https://github.com/sergiocorreia/reghdfe).

The command handles panel models applied in applied econometrics:
- **Class A**: Individual + time FE (e.g., firm panel)
- **Class B**: Exporter-year + importer-year FE (two-way gravity)
- **Class C**: Exporter-year + importer-year + pair FE (three-way gravity)

## Installation

### From GitHub (recommended)

```stata
net install ivppmlhdfe, from("https://raw.githubusercontent.com/ekwonomist/ivppmlhdfe/main/") replace
```

### Dependencies

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

```stata
* Two-way gravity model
ivppmlhdfe trade (policy = instrument), absorb(exp_year imp_year) vce(cluster pair_id)

* With bias correction (Class A)
ivppmlhdfe y (x = z), absorb(id year) vce(cluster id) biascorrection(a id year)

* Three-way gravity model with pair FE
ivppmlhdfe trade (policy = instrument), absorb(exp_year imp_year pair_id) vce(cluster pair_id)

* Incidence-rate ratios
ivppmlhdfe trade (policy = instrument), absorb(exp_year imp_year) vce(cluster pair_id) irr
```

## References

- Mullahy, J. (1997). "Instrumental-variable estimation of count data models." *Review of Economics and Statistics*, 79(4), 586-593.
- Correia, S., P. Guimaraes, and T. Zylkin (2020). "Fast Poisson estimation with high-dimensional fixed effects." *Stata Journal*, 20(1), 95-115.
- Weidner, M. and T. Zylkin (2021). "Bias and consistency in three-way gravity models." *Journal of International Economics*, 132, 103513.
- Fernandez-Val, I. and M. Weidner (2016). "Individual and time effects in nonlinear panel models with large N, T." *Journal of Econometrics*, 192(1), 291-312.