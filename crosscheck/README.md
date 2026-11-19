# Cross-language check

The four versions of `ivppmlhdfe` (Stata/Mata, Julia, R, Python) are tested on one shared set of 17
cases. This folder holds the data, the case list, every version's results, and the script that
compares them, so the agreement can be checked again from a clone of the repository.

## Contents

| Path | What it is |
|---|---|
| `testdata/cases.json` | The 17 cases: depvar, regressors, instruments, `absorb()`, VCE, weights, offset |
| `testdata/class*.csv` | Simulated panels (Class A, B, C, and Class A with singletons), 17 significant digits |
| `testdata/make_cases.py` | Rebuilds the CSVs and `cases.json` (numpy, fixed seed) |
| `results/<lang>_accel.json` | Results with the fast path on (the default) |
| `results/<lang>_noaccel.json` | Results with the fast path off |
| `results/stata_v094.json` | Stata v0.9.4, the reference for the earlier algorithm |
| `results/julia_v094.json` | Julia v0.9.4 |
| `results/stata_ppmlhdfe.json` | `ppmlhdfe` on the two cases without instruments |
| `results/stata_ffirst_reference.json` | First-stage statistics from `ivreghdfe`/`ranktest` at the final IRLS weights |
| `SPEC_ffirst.md` | The first-stage formulas every version implements |
| `compare.py` | Compares two result files |

The cases cover Class A (`id year`), B (`exp#year imp#year`) and C (adds `exp#imp`); robust, one-way
and two-way clustered SEs; just-identified, overidentified and two-endogenous models; pweights and
fweights; an offset; PPML without instruments; separation (all-zero groups); and singletons.

## Comparing two versions

```
python crosscheck/compare.py crosscheck/results/r_accel.json crosscheck/results/stata_accel.json --ffirst
```

The difference is |a - b| / max(|b|, 1e-3), reported per case for coefficients, SEs and (with
`--ffirst`) first-stage statistics. N, singletons and separated observations must match exactly.
The script exits with code 1 when a case is outside `--tol` (default 1e-6).

## Results at v1.0.0 (largest difference over the 17 cases)

| Comparison | Coefficients | SEs | First stage |
|---|---|---|---|
| Stata fast path off vs Stata v0.9.4 | 0 (bit-identical) | 0 | |
| Julia fast path off vs Julia v0.9.4 | 0 (bit-identical, 16 of 17 cases; see note) | 0 | |
| Stata fast path on vs off | 3.2e-11 | 6.2e-9 | 8.5e-9 |
| R, Python vs Stata, fast path on | 1.6e-15 | 1.6e-15 | 7.4e-15 |
| Julia vs Stata, fast path on | 7.2e-14 | 4.4e-10 | 2.3e-10 |
| R, Python fast path off vs Stata v0.9.4 | 5.1e-11 | 1.7e-8 | |
| Julia fast path off vs Stata v0.9.4 | 1.4e-10 | 6.5e-7 | |
| Each version's first stage vs `ivreghdfe` reference | | | 1.5e-9 |

Note: Julia v0.9.4 treated frequency weights as probability weights. In v1.0.0 the keyword
`weight_type=:fweight` gives Stata's frequency-weight results; `A_fw` therefore differs from
`julia_v094.json` (N = 1999 against 1000) and matches Stata.

## Running the package tests on these cases

- Python: `pytest python/tests` finds this folder by itself.
- R: set `IVPPML_TESTDATA` to `crosscheck/testdata`, `IVPPML_REFERENCE` to
  `crosscheck/results/stata_v094.json` and `IVPPML_FFIRST_REFERENCE` to
  `crosscheck/results/stata_ffirst_reference.json`, then run `testthat::test_local("R/ivppmlhdfe")`.
  Without them the case tests are skipped and the package's own unit tests run.
