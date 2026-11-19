"""Build the shared cross-language test cases for ivppmlhdfe v1.0.0.

Reads the three example datasets in github_repo/data, adds deterministic extra columns
(second and third instruments, a second endogenous regressor, weights, an offset, zeros
with a separated group, singletons), writes one CSV per dataset with 17 significant
digits (doubles round-trip exactly), and writes cases.json.

Run from the project root:  py -3.12 program/sandbox/v010_build/testdata/make_cases.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
SRC = ROOT / "github_repo" / "data"
OUT = Path(__file__).resolve().parent


def augment(d: pd.DataFrame, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    n = len(d)
    d = d.copy()
    d["z2"] = rng.standard_normal(n)
    d["z3"] = rng.standard_normal(n)
    # second endogenous regressor: driven by z2, z3 and by x1's first-stage error
    d["x3"] = 0.5 * d["z2"] + 0.3 * d["z3"] + 0.4 * (d["x1"] - 0.8 * d["z"]) + 0.5 * rng.standard_normal(n)
    d["w"] = rng.uniform(0.5, 2.0, n)              # pweight
    d["fw"] = 1 + (np.arange(n) % 3)               # fweight 1..3
    d["off"] = 0.2 * rng.standard_normal(n)        # offset
    y0 = d["y"].to_numpy().copy()
    y0[rng.uniform(size=n) < 0.25] = 0.0           # 25% zeros
    d["y0"] = y0
    return d


def write(d: pd.DataFrame, name: str) -> str:
    f = OUT / f"{name}.csv"
    d.to_csv(f, index=False, float_format="%.17g", lineterminator="\n")
    return f.name


def main() -> None:
    A = augment(pd.read_stata(SRC / "ivppmlhdfe_ClassA.dta"), 101)
    B = augment(pd.read_stata(SRC / "ivppmlhdfe_ClassB.dta"), 202)
    C = augment(pd.read_stata(SRC / "ivppmlhdfe_ClassC.dta"), 303)

    # B with a separated exporter-year group: all y0 = 0 for exp 3 in year 5
    B.loc[(B["exp"] == 3) & (B["year"] == 5), "y0"] = 0.0

    # A with 5 singleton observations (new ids seen once)
    extra = A.sample(5, random_state=7).copy()
    extra["id"] = np.arange(1001, 1006)
    A_single = pd.concat([A, extra], ignore_index=True)

    files = {"A": write(A, "classA"), "B": write(B, "classB"), "C": write(C, "classC"),
             "A_single": write(A_single, "classA_single")}

    fe_A = [["id"], ["year"]]
    fe_B = [["exp", "year"], ["imp", "year"]]
    fe_C = [["exp", "imp"], ["exp", "year"], ["imp", "year"]]
    rob = {"type": "robust"}

    def cl(*v):
        return {"type": "cluster", "vars": list(v)}

    def case(cid, data, absorb, vce, endog=("x1",), instruments=("z",), exog=("x2",),
             depvar="y", weights=None, weight_type=None, offset=None, note=""):
        return {"id": cid, "data": files[data], "depvar": depvar, "exog": list(exog),
                "endog": list(endog), "instruments": list(instruments), "absorb": absorb,
                "vce": vce, "weights": weights, "weight_type": weight_type,
                "offset": offset, "note": note}

    cases = [
        case("A_rob", "A", fe_A, rob),
        case("A_cl", "A", fe_A, cl("id")),
        case("B_rob", "B", fe_B, rob),
        case("B_cl", "B", fe_B, cl("pair")),
        case("C_rob", "C", fe_C, rob),
        case("C_cl", "C", fe_C, cl("pair")),
        case("C_cl2", "C", fe_C, cl("exp", "imp"), note="two-way clustering"),
        case("A_overid", "A", fe_A, rob, instruments=("z", "z2"), note="overidentified"),
        case("C_overid", "C", fe_C, cl("pair"), instruments=("z", "z2"), note="overidentified"),
        case("C_k2", "C", fe_C, cl("pair"), endog=("x1", "x3"), instruments=("z", "z2", "z3"),
             note="two endogenous regressors, overidentified"),
        case("B_pw", "B", fe_B, cl("pair"), weights="w", weight_type="pweight"),
        case("A_fw", "A", fe_A, rob, weights="fw", weight_type="fweight"),
        case("C_off", "C", fe_C, rob, offset="off"),
        case("A_ppml", "A", fe_A, rob, endog=(), instruments=(), exog=("x1", "x2"),
             note="no instruments: must match ppmlhdfe"),
        case("C_ppml", "C", fe_C, cl("pair"), endog=(), instruments=(), exog=("x1", "x2"),
             note="no instruments: must match ppmlhdfe"),
        case("B_sep", "B", fe_B, rob, depvar="y0",
             note="25% zeros; exporter 3 in year 5 is all zero, so separated by the fe method"),
        case("A_single", "A_single", fe_A, rob, note="5 singleton observations to drop"),
    ]
    (OUT / "cases.json").write_text(json.dumps({"version": "1.0.0", "cases": cases}, indent=2) + "\n",
                                    encoding="utf-8", newline="\n")
    for k, v in files.items():
        print(k, v, (OUT / v).stat().st_size)
    print(len(cases), "cases written to", OUT / "cases.json")


if __name__ == "__main__":
    main()
