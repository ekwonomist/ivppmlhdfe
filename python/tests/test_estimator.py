"""Estimator checks: analytic cases, accel on vs off, separation, options."""
import json
import warnings

import numpy as np
import pandas as pd
import pytest

import ivppmlhdfe
from ivppmlhdfe import FESolver
from ivppmlhdfe._fit import _tsls

from conftest import RESULTS, case_args


def _rel(a, b):
    return abs(a - b) / max(abs(b), 1e-3)


# ---------------------------------------------------------------------------------------------
# analytic checks
# ---------------------------------------------------------------------------------------------
def _small_panel(seed=0, n=12, T=8):
    rng = np.random.default_rng(seed)
    i = np.repeat(np.arange(n), T)
    t = np.tile(np.arange(T), n)
    a = 0.3 * rng.standard_normal(n)[i] + 0.2 * rng.standard_normal(T)[t]
    z = rng.standard_normal(n * T)
    x2 = 0.5 * rng.standard_normal(n * T)
    e = rng.standard_normal(n * T)
    x1 = 0.8 * z + e + 0.4 * a
    y = np.exp(0.5 * x1 + 0.3 * x2 + a + 0.5 * e * 0.5 + 0.3 * rng.standard_normal(n * T))
    return pd.DataFrame({"id": i, "year": t, "z": z, "x1": x1, "x2": x2, "y": y})


def _newton_poisson(y, X):
    b = np.zeros(X.shape[1])
    b[0] = np.log(y.mean())
    for _ in range(200):
        mu = np.exp(X @ b)
        H = X.T @ (mu[:, None] * X)
        step = np.linalg.solve(H, X.T @ (y - mu))
        b += step
        if np.max(np.abs(step)) < 1e-13:
            break
    mu = np.exp(X @ b)
    H = X.T @ (mu[:, None] * X)
    s = X * (y - mu)[:, None]
    Hi = np.linalg.inv(H)
    N = len(y)
    V = N / (N - 1) * Hi @ (s.T @ s) @ Hi
    return b, V


def test_ppml_matches_newton_on_dummies():
    df = _small_panel(1)
    r = ivppmlhdfe.fit(df, "y ~ x1 + x2 | id + year", tol=1e-12)
    Did = np.eye(12)[df["id"]]
    Dyr = np.eye(8)[df["year"]][:, 1:]
    X = np.column_stack([Did, Dyr, df["x1"], df["x2"]])
    b, V = _newton_poisson(df["y"].to_numpy(), X)
    assert abs(r.coef["x1"] - b[-2]) < 1e-9 and abs(r.coef["x2"] - b[-1]) < 1e-9
    se = np.sqrt(np.diag(V))
    assert _rel(r.se["x1"], se[-2]) < 1e-7 and _rel(r.se["x2"], se[-1]) < 1e-7


def test_iv_step_is_2sls_with_dummies_at_fixed_weights():
    rng = np.random.default_rng(2)
    N = 360
    c1 = np.repeat(np.arange(30), 12)
    c2 = np.tile(np.arange(12), 30)
    c3 = rng.integers(0, 6, N)
    c3[:6] = np.arange(6)
    w = rng.uniform(0.2, 3.0, N)
    zdep = rng.standard_normal(N)
    X = rng.standard_normal((N, 2))          # [exog x2, endog x1]
    inst = rng.standard_normal((N, 1))
    S = FESolver([c1, c2, c3], [30, 12, 6]).update_weights(w)
    M, _, conv = S.partial_out(np.column_stack([zdep, X, inst]), tol=1e-14)
    assert conv.all()
    b, _, _ = _tsls(M[:, 0], M[:, 1:3], np.column_stack([M[:, 1], M[:, 3]]), w)
    D = np.concatenate([np.eye(30)[c1], np.eye(12)[c2][:, 1:], np.eye(6)[c3][:, 1:]], axis=1)
    XX = np.column_stack([X, D])
    ZZ = np.column_stack([X[:, :1], inst, D])
    Zw = ZZ * w[:, None]
    Pi = np.linalg.solve(Zw.T @ ZZ, Zw.T @ XX)
    Xh = ZZ @ Pi
    bd = np.linalg.solve((Xh * w[:, None]).T @ XX, (Xh * w[:, None]).T @ zdep)
    assert np.max(np.abs(b - bd[:2])) < 1e-9


def test_iv_moment_conditions_hold_at_convergence(cases, datasets):
    c = next(c for c in cases if c["id"] == "C_rob")
    df = datasets(c["data"])
    f, kw = case_args(c)
    r = ivppmlhdfe.fit(df, f, tol=1e-12, **kw)
    mu = r.predict("mu").to_numpy()
    y = df["y"].to_numpy()
    u = y - mu
    for q in ("x2", "z"):
        v = df[q].to_numpy()
        assert abs(np.sum(v * u)) / np.sum(np.abs(v * u)) < 1e-8
    for fe in (["exp", "imp"], ["exp", "year"], ["imp", "year"]):
        g = df.groupby(fe).ngroup().to_numpy()
        s = np.bincount(g, weights=u)
        assert np.max(np.abs(s)) / np.sum(np.abs(u)) < 1e-8


# ---------------------------------------------------------------------------------------------
# the shared cases
# ---------------------------------------------------------------------------------------------
def test_all_cases_accel_on_vs_off(cases, datasets):
    bad = []
    for c in cases:
        f, kw = case_args(c)
        df = datasets(c["data"])
        a = ivppmlhdfe.fit(df, f, accel=True, **kw)
        b = ivppmlhdfe.fit(df, f, accel=False, **kw)
        assert a.converged and b.converged
        assert (a.N, a.num_singletons, a.num_separated) == (b.N, b.num_singletons, b.num_separated)
        rc = max(_rel(a.coef[k], b.coef[k]) for k in a.coef.index)
        rs = max(_rel(a.se[k], b.se[k]) for k in a.se.index)
        if rc > 1e-8 or rs > 1e-6:
            bad.append((c["id"], rc, rs))
    assert not bad, bad


def test_all_cases_vs_stata_v094(cases, datasets):
    ref = RESULTS / "stata_v094.json"
    if not ref.exists():
        pytest.skip("results/stata_v094.json not available")
    S = json.loads(ref.read_text(encoding="utf-8"))["cases"]
    bad = []
    for c in cases:
        if c["id"] not in S or not S[c["id"]].get("ok", True):
            continue
        s = S[c["id"]]
        f, kw = case_args(c)
        for accel in (True, False):
            r = ivppmlhdfe.fit(datasets(c["data"]), f, accel=accel, **kw)
            counts = (r.N, r.num_singletons, r.num_separated) == \
                (s["N"], s["num_singletons"], s["num_separated"])
            rc = max(_rel(r.coef[k], s["coef"][k]) for k in s["coef"])
            rs = max(_rel(r.se[k], s["se"][k]) for k in s["se"])
            if not counts or rc > 1e-6 or rs > 1e-6:
                bad.append((c["id"], accel, counts, rc, rs))
    assert not bad, bad


def _flat(ff):
    out = {}
    for k, v in (ff or {}).items():
        if isinstance(v, dict):
            out.update({f"{k}.{s}": x for s, x in v.items()})
        else:
            out[k] = v
    return out


def test_ffirst_vs_stata_reference(cases, datasets):
    ref = RESULTS / "stata_ffirst_reference.json"
    if not ref.exists():
        pytest.skip("results/stata_ffirst_reference.json not available")
    S = json.loads(ref.read_text(encoding="utf-8"))["cases"]
    bad = []
    for c in cases:
        if c["id"] not in S or not c["endog"]:
            continue
        f, kw = case_args(c)
        r = ivppmlhdfe.fit(datasets(c["data"]), f, accel=True, ffirst=True, **kw)
        a, b = _flat(r.ffirst), _flat(S[c["id"]]["ffirst"])
        assert set(a) == set(b), (c["id"], set(a) ^ set(b))
        worst = max(_rel(a[k], b[k]) for k in b)
        if worst > 1e-6:
            bad.append((c["id"], worst))
    assert not bad, bad


# ---------------------------------------------------------------------------------------------
# separation, singletons, options
# ---------------------------------------------------------------------------------------------
def test_counts_singletons_and_fe_separation(cases, datasets):
    by = {c["id"]: c for c in cases}
    f, kw = case_args(by["A_single"])
    r = ivppmlhdfe.fit(datasets("classA_single.csv"), f, **kw)
    assert r.num_singletons == 5 and r.N == 1000
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rk = ivppmlhdfe.fit(datasets("classA_single.csv"), f, keep_singletons=True, **kw)
    assert rk.num_singletons == 0 and rk.N == 1005
    f, kw = case_args(by["B_sep"])
    r = ivppmlhdfe.fit(datasets("classB.csv"), f, **kw)
    # Stata convention: FE-separated obs are counted with the singletons (reghdfe iweight trick)
    assert r.num_singletons == 29 and r.num_separated == 0 and r.N == 8700 - 29


def test_relu_detects_regressor_separation():
    df = _small_panel(3, n=20, T=10)
    rng = np.random.default_rng(9)
    d = np.zeros(len(df))
    d[rng.choice(len(df), 7, replace=False)] = 1.0
    df["d"] = d
    df.loc[df["d"] == 1, "y"] = 0.0
    df.loc[rng.choice(np.flatnonzero(d == 0), 15, replace=False), "y"] = 0.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = ivppmlhdfe.fit(df, "y ~ x2 + d | id + year | x1 ~ z")
        r0 = ivppmlhdfe.fit(df[df["d"] == 0], "y ~ x2 | id + year | x1 ~ z")
    assert r.num_separated == 7
    assert "d" not in r.coef.index
    assert abs(r.coef["x1"] - r0.coef["x1"]) < 1e-9


def test_formula_and_errors():
    p = ivppmlhdfe.parse_formula("y ~ x2 + x4 | id + exp^year | x1 + x3 ~ z + z2")
    assert p.exog == ["x2", "x4"] and p.endog == ["x1", "x3"] and p.instruments == ["z", "z2"]
    assert p.absorb == [["id"], ["exp", "year"]]
    p = ivppmlhdfe.parse_formula("y ~ x2 | x1 ~ z")
    assert p.absorb == [] and p.endog == ["x1"]
    df = _small_panel(4)
    with pytest.raises(ValueError):
        ivppmlhdfe.fit(df.assign(y=-df["y"]), "y ~ x2 | id")
    with pytest.raises(ValueError):
        ivppmlhdfe.fit(df, "y ~ 1 | id")
    with pytest.raises(ValueError):
        ivppmlhdfe.fit(df, "y ~ x2 | id | x1 ~ x1")


def test_noabsorb_and_predict():
    df = _small_panel(5)
    r = ivppmlhdfe.fit(df, "y ~ x2 | x1 ~ z")
    mu = r.predict("mu")
    assert abs(mu.sum() - df["y"].sum()) / df["y"].sum() < 1e-8   # intercept FOC
    assert np.allclose(np.log(mu), r.predict("xb") + r.predict("d"))


def test_collinear_regressor_keeps_the_earlier_one():
    df = _small_panel(7).assign(x2b=lambda d: 2.0 * d["x2"])
    with pytest.warns(RuntimeWarning, match="x2b"):
        r = ivppmlhdfe.fit(df, "y ~ x2 + x2b | id + year | x1 ~ z")
    r0 = ivppmlhdfe.fit(df, "y ~ x2 | id + year | x1 ~ z")
    assert list(r.coef.index) == ["x2", "x1"]
    assert np.allclose(r.coef.to_numpy(), r0.coef.to_numpy(), rtol=1e-12, atol=0)


def test_cluster_single_group_demotes_to_robust():
    df = _small_panel(6).assign(g=1)
    with pytest.warns(RuntimeWarning):
        r = ivppmlhdfe.fit(df, "y ~ x2 | id + year | x1 ~ z", vcov={"CL": ["g"]})
    assert r.vce_type == "robust"
