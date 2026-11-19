"""ivppmlhdfe.fit -- IV-PPML with high-dimensional fixed effects (IRLS-IV).

The additive moment E[q (y - mu)] = 0 with q = (x', z')' and mu = exp(x'b + sum of FE + offset),
solved by iteratively reweighted 2SLS with the fixed effects concentrated out at each iteration
(Windmeijer and Santos Silva 1997; Correia, Guimaraes and Zylkin 2020). Everything that is not
about the fixed-effect solver reproduces Stata ivppmlhdfe v0.9.4 (see the package README and
REPORT.md for the line-by-line map).

Fast path (``accel=True``, the default; SPEC section 2):
  W  warm start: partial out (z~_{k-1} + z_k - z_{k-1}, X~_{k-1}, Z~_{k-1}); iteration 1 starts X, Z
     from the collinearity pass;
  T  fereg's inner-tolerance schedule (8 log-linear stages 1e-4 -> target, 0.1*eps, the target once
     eps < tol, never loosened);
  R  VCE refresh: after convergence X, Z (and z) are partialled out once more at the converged
     weights, and e(V) and ffirst are computed from those columns.
The FE solver (S: exact elimination of the largest FE + Jacobi-PCG) is used in both modes.
``accel=False``: cold start, ppmlhdfe's tolerance rule, v0.9.4's VCE (last iteration's columns).
"""
from __future__ import annotations

import time
import warnings

import numpy as np
import pandas as pd
from scipy.special import gammaln

from ._fe import FESolver
from ._formula import parse_formula
from ._prep import (compact_codes, drop_fe_stage, fe_codes, relu_separation,
                    select_not_collinear)
from ._result import IVPPMLResult
from ._vce import sandwich_vce

__all__ = ["fit", "IVPPMLError", "IVPPMLConvergenceError"]

VERSION = "1.0.0"
EPS100 = float(np.spacing(100.0))      # 1.42e-14, ppmlhdfe's floor for mu


class IVPPMLError(RuntimeError):
    """Estimation failure (Stata rc 9003: runaway divergence)."""


class IVPPMLConvergenceError(IVPPMLError):
    """IRLS did not converge within maxit (Stata rc 430)."""


# ----------------------------------------------------------------------------------------------
# option parsing
# ----------------------------------------------------------------------------------------------
def _parse_vcov(vcov):
    if vcov is None:
        return "robust", []
    if isinstance(vcov, str):
        v = vcov.strip()
        low = v.lower()
        if low in ("robust", "hetero", "hc1", "hc0", "r", "unadjusted", "ols", "conventional"):
            return "robust", []
        for pre in ("cl:", "cluster:", "cluster ", "cl "):
            if low.startswith(pre):
                vars_ = [t for t in v[len(pre):].replace("+", " ").split() if t]
                return "cluster", vars_
        raise ValueError(f"ivppmlhdfe: unknown vcov {vcov!r}; use 'robust' or {{'CL': ['var']}}")
    if isinstance(vcov, dict):
        if len(vcov) != 1:
            raise ValueError("ivppmlhdfe: vcov dict must have one key, 'CL'")
        (k, v), = vcov.items()
        if str(k).upper() not in ("CL", "CLUSTER"):
            raise ValueError("ivppmlhdfe: vcov dict key must be 'CL'")
        if isinstance(v, str):
            v = [t for t in v.replace("+", " ").split() if t]
        v = list(v)
        if not v:
            raise ValueError("ivppmlhdfe: vcov={'CL': [...]} needs at least one variable")
        return "cluster", v
    if isinstance(vcov, (list, tuple)):
        return "cluster", list(vcov)
    raise TypeError("ivppmlhdfe: vcov must be 'robust' or {'CL': [vars]}")


def _parse_separation(sep):
    if sep is None or sep is False:
        return set()
    if isinstance(sep, str):
        s = sep.strip().lower()
        if s in ("", "def", "default", "standard", "on", "auto"):
            return {"fe", "relu"}
        if s in ("all", "full"):
            return {"fe", "relu", "mu"}
        if s in ("no", "off", "none"):
            return set()
        sep = s.replace(",", " ").split()
    out = set()
    for t in sep:
        t = str(t).strip().lower()
        if t == "ir":
            t = "relu"
        if t not in ("fe", "simplex", "relu", "mu"):
            raise ValueError(f"ivppmlhdfe: unknown separation method {t!r} (fe, simplex, relu, mu)")
        out.add(t)
    return out


def _parse_accel(accel):
    """True -> W+T+R, False -> none; a string of letters from 'wtr' picks components."""
    if accel is True:
        return {"warmstart": True, "tolpolicy": "fereg", "refresh": True}
    if accel is False or accel is None:
        return {"warmstart": False, "tolpolicy": "ppmlhdfe", "refresh": False}
    if isinstance(accel, str):
        s = accel.lower()
        bad = set(s) - set("wtrs ,")
        if bad:
            raise ValueError("ivppmlhdfe: accel string may contain only w, t, r (and s)")
        return {"warmstart": "w" in s, "tolpolicy": "fereg" if "t" in s else "ppmlhdfe",
                "refresh": "r" in s}
    if isinstance(accel, dict):
        d = _parse_accel(False)
        d.update(accel)
        return d
    raise TypeError("ivppmlhdfe: accel must be True, False or a string such as 'wt'")


def _num_col(df, name, role):
    if name not in df.columns:
        raise KeyError(f"ivppmlhdfe: {role} column {name!r} not found")
    s = df[name]
    try:
        a = pd.to_numeric(s, errors="raise").to_numpy(dtype=np.float64, na_value=np.nan)
    except (TypeError, ValueError) as e:
        raise TypeError(f"ivppmlhdfe: {role} column {name!r} is not numeric") from e
    if np.isinf(a).any():
        raise ValueError(f"ivppmlhdfe: {role} column {name!r} contains infinite values")
    return a


# ----------------------------------------------------------------------------------------------
# 2SLS
# ----------------------------------------------------------------------------------------------
def _fstack(*blocks):
    """Column-stack 1-d / 2-d blocks into a Fortran-ordered (column-contiguous) array."""
    blocks = [b[:, None] if b.ndim == 1 else b for b in blocks]
    N = blocks[0].shape[0]
    out = np.empty((N, sum(b.shape[1] for b in blocks)), order="F")
    j = 0
    for b in blocks:
        out[:, j:j + b.shape[1]] = b
        j += b.shape[1]
    return out


def _solve_sym(A, B):
    try:
        return np.linalg.solve(A, B)
    except np.linalg.LinAlgError:
        return np.linalg.pinv(A) @ B


def _tsls(z_dm, X_dm, Z_dm, w):
    """Weighted 2SLS without constant: b, resid, Xhat."""
    Zw = Z_dm * w[:, None]
    ZwZ = Zw.T @ Z_dm
    ZwX = Zw.T @ X_dm
    Pi = _solve_sym(ZwZ, ZwX)
    Xhat = Z_dm @ Pi
    XhwX = Pi.T @ ZwX
    Xhwy = Pi.T @ (Zw.T @ z_dm)
    b = _solve_sym(0.5 * (XhwX + XhwX.T), Xhwy)
    resid = z_dm - X_dm @ b
    return b, resid, Xhat


# ----------------------------------------------------------------------------------------------
# main entry point
# ----------------------------------------------------------------------------------------------
def fit(data, formula, vcov="robust", weights=None, weight_type="pweight", offset=None,
        exposure=None, accel=True, ffirst=False, tol=1e-8, itol=None, maxit=10000,
        separation=("fe", "relu"), keep_singletons=False, verbose=0, guess="simple"):
    """Estimate an IV-PPML model with high-dimensional fixed effects.

    Parameters
    ----------
    data : pandas.DataFrame
    formula : str
        ``"y ~ x2 | id + year | x1 ~ z"`` (fixest style; interactions ``exp^year``; the FE part
        and the IV part are optional).
    vcov : "robust" or {"CL": ["var", ...]}
        Robust (default) or one-, two- or three-way clustered standard errors.
    weights : str or None
        Name of the weight column; ``weight_type`` is "pweight" (default), "aweight" or "fweight".
    offset, exposure : str or None
        Offset column, or exposure column (offset = log(exposure)).
    accel : bool or str
        Fast path on (default) or off; a string of 'w' (warm start), 't' (fereg tolerance
        schedule), 'r' (VCE refresh) picks components (diagnostics).
    ffirst : bool
        First-stage statistics (SPEC section 3).
    tol, itol, maxit : convergence tolerance (deviance), inner FE tolerance target
        (default max(1e-12, 0.1*tol)), maximum IRLS iterations.
    separation : tuple of "fe", "relu", "simplex", "mu" (or "default", "all", "none")
    keep_singletons : bool
    verbose : int
        0 silent, 1 iteration log, 2 more.
    guess : "simple" or "mean"
    """
    t_fit0 = time.perf_counter()
    flags = _parse_accel(accel)
    pf = parse_formula(formula)
    vce_type, clustvars = _parse_vcov(vcov)
    sep = _parse_separation(separation)
    do_sep_fe = "fe" in sep
    do_sep_relu = ("relu" in sep) or ("simplex" in sep)
    do_sep_mu = "mu" in sep
    if weight_type not in ("pweight", "aweight", "fweight"):
        raise ValueError("ivppmlhdfe: weight_type must be 'pweight', 'aweight' or 'fweight'")
    if offset is not None and exposure is not None:
        raise ValueError("ivppmlhdfe: cannot specify both offset and exposure")
    if guess not in ("simple", "mean", "default"):
        raise ValueError("ivppmlhdfe: guess must be 'simple' or 'mean'")
    if not isinstance(data, pd.DataFrame):
        data = pd.DataFrame(data)
    index = data.index
    N_rows_data = len(data)

    exog, endog, inst = list(pf.exog), list(pf.endog), list(pf.instruments)
    if len(exog) + len(endog) == 0:
        raise ValueError("ivppmlhdfe: specify at least one regressor (an empty model is not supported)")
    if len(inst) < len(endog):
        raise ValueError("ivppmlhdfe: equation not identified (fewer instruments than endogenous regressors)")
    bad = set(inst) & set(endog)
    if bad:
        raise ValueError(f"ivppmlhdfe: instrument cannot equal endogenous regressor: {sorted(bad)}")
    bad = set(inst) & set(exog)
    if bad:
        raise ValueError(f"ivppmlhdfe: instrument cannot equal exogenous regressor: {sorted(bad)}")

    # ---- markout ---------------------------------------------------------------------------
    y_all = _num_col(data, pf.depvar, "dependent variable")
    touse = ~np.isnan(y_all)
    num = {}
    for c in exog + endog + inst:
        num[c] = _num_col(data, c, "regressor/instrument")
        touse &= ~np.isnan(num[c])
    w_all = None
    if weights is not None:
        w_all = _num_col(data, weights, "weight")
        touse &= ~np.isnan(w_all)
    off_all = None
    if exposure is not None:
        ex = _num_col(data, exposure, "exposure")
        touse &= ~np.isnan(ex)
        if np.any(ex[touse] <= 0):
            raise ValueError("ivppmlhdfe: exposure must be greater than zero")
        with np.errstate(divide="ignore"):
            off_all = np.log(ex)
    elif offset is not None:
        off_all = _num_col(data, offset, "offset")
        touse &= ~np.isnan(off_all)
    fe_cols = sorted({c for t in pf.absorb for c in t})
    cl_terms = [[p.strip() for p in c.replace("#", "^").split("^")] for c in clustvars]
    for c in fe_cols + sorted({p for t in cl_terms for p in t}):
        if c not in data.columns:
            raise KeyError(f"ivppmlhdfe: column {c!r} not found")
        touse &= ~pd.isna(data[c]).to_numpy()
    if np.any(y_all[touse] < 0):
        raise ValueError(f"ivppmlhdfe: {pf.depvar} must be nonnegative")
    if w_all is not None:
        if np.any(w_all[touse] < 0):
            raise ValueError("ivppmlhdfe: negative weights encountered")
        touse &= ~(w_all == 0)
    rows = np.flatnonzero(touse)
    N_full = int(rows.size)
    if N_full == 0:
        raise ValueError("ivppmlhdfe: no observations")
    for c in inst:
        v = num[c][rows]
        if v.min() == v.max():
            raise ValueError(f"ivppmlhdfe: instrument {c} has zero variance (constant column)")

    # ---- cluster codes (G = 1 demotion before the FE stage, as in Stata) --------------------
    cl_codes = []
    if vce_type == "cluster":
        cl_codes, cl_nlev = fe_codes(data, cl_terms, rows)
        if len(cl_codes) == 1 and cl_nlev[0] <= 1:
            warnings.warn(f"ivppmlhdfe: only 1 cluster in {clustvars[0]}; switching to vce(robust)",
                          RuntimeWarning, stacklevel=2)
            vce_type, clustvars, cl_codes, cl_terms = "robust", [], [], []

    # ---- fixed-effect codes and the FE stage (singletons + FE separation) ------------------
    absorb = [list(t) for t in pf.absorb]
    if absorb:
        codes, nlev = fe_codes(data, absorb, rows)
    else:
        codes, nlev = [np.zeros(N_full, dtype=np.int64)], [1]       # noabsorb: one constant FE
    y = y_all[rows]
    if keep_singletons:
        warnings.warn("ivppmlhdfe: keeping singleton groups will keep fixed effects that cause "
                      "separation", RuntimeWarning, stacklevel=2)
        keep = np.ones(N_full, dtype=bool)
    else:
        keep = drop_fe_stage(codes, nlev, np.ones(N_full, dtype=bool),
                             y=y if do_sep_fe else None)
    num_singletons = int(N_full - keep.sum())
    if keep.sum() == 0:
        raise ValueError("ivppmlhdfe: no observations remaining after dropping singletons and "
                         "separated observations")

    def _subset(mask):
        nonlocal rows, codes, nlev, cl_codes
        rows = rows[mask]
        codes, nlev = compact_codes(codes, mask)
        cl_codes = [c[mask] for c in cl_codes]

    if not keep.all():
        _subset(keep)

    def _load():
        y_ = y_all[rows]
        X_ = _fstack(*[num[c][rows] for c in exog + endog]) if exog + endog else \
            np.zeros((rows.size, 0), order="F")
        Z_ = _fstack(*[num[c][rows] for c in inst]) if inst else np.zeros((rows.size, 0), order="F")
        w_ = w_all[rows] if w_all is not None else np.ones(rows.size)
        o_ = off_all[rows] if off_all is not None else np.zeros(rows.size)
        return y_, X_, Z_, w_, o_

    y, X, Z, w_user, off_orig = _load()
    if y.sum() == 0:
        raise ValueError(f"ivppmlhdfe: {pf.depvar} is identically zero in the estimation sample")

    tol = float(tol)
    target_inner_tol = float(itol) if (itol is not None and itol > 0) else max(1e-12, 0.1 * tol)
    stats = {"t_coll": 0.0, "t_setup": 0.0, "t_partial": 0.0, "t_relu": 0.0, "t_refresh": 0.0,
             "inner_iters": 0, "passes": 0, "solves": 0, "inner_allconv": True,
             "tol_trace": [], "it_trace": []}

    # ---- collinearity (partial out under the user weights at the target tolerance) ----------
    def _collinearity(solver):
        nonlocal exog, endog, inst, X, Z
        t0 = time.perf_counter()
        solver.update_weights(w_user)
        raw = _fstack(X, Z)
        part, its, cv = solver.partial_out(raw, target_inner_tol) if raw.shape[1] else (raw, [], [])
        stats["inner_iters"] += int(np.sum(its))
        stats["inner_allconv"] &= bool(np.all(cv))
        Kx, Kd, L0 = len(exog), len(endog), len(inst)
        ke = select_not_collinear(part[:, :Kx], raw[:, :Kx], w_user)
        kd = select_not_collinear(part[:, Kx:Kx + Kd], raw[:, Kx:Kx + Kd], w_user)
        ki = select_not_collinear(part[:, Kx + Kd:], raw[:, Kx + Kd:], w_user)
        dropped = (Kx - len(ke), Kd - len(kd), L0 - len(ki))
        cols = list(ke) + [Kx + j for j in kd] + [Kx + Kd + j for j in ki]
        if sum(dropped):
            names = [exog[j] for j in range(Kx) if j not in ke] + \
                    [endog[j] for j in range(Kd) if j not in kd] + \
                    [inst[j] for j in range(L0) if j not in ki]
            warnings.warn(f"ivppmlhdfe: omitted because of collinearity: {' '.join(names)}",
                          RuntimeWarning, stacklevel=3)
            exog = [exog[j] for j in ke]
            endog = [endog[j] for j in kd]
            inst = [inst[j] for j in ki]
            X = np.asfortranarray(raw[:, list(ke) + [Kx + j for j in kd]])
            Z = np.asfortranarray(raw[:, [Kx + Kd + j for j in ki]])
            if len(inst) < len(endog):
                raise ValueError("ivppmlhdfe: equation not identified after collinearity removal")
            if len(exog) + len(endog) == 0:
                raise ValueError("ivppmlhdfe: no regressor left after collinearity removal")
        stats["t_coll"] += time.perf_counter() - t0
        return np.asfortranarray(part[:, cols])

    t0 = time.perf_counter()
    solver = FESolver(codes, nlev)
    stats["t_setup"] += time.perf_counter() - t0
    coll_part = _collinearity(solver)

    # ---- ReLU separation ----------------------------------------------------------------------
    num_sep_relu = 0
    if do_sep_relu and X.shape[1] > 0 and np.any(y == 0):
        t0 = time.perf_counter()
        sep_idx = relu_separation(solver, y, X, target_inner_tol, verbose=verbose)
        if sep_idx.size:
            m = np.ones(rows.size, dtype=bool)
            m[sep_idx] = False
            m = drop_fe_stage(codes, nlev, m,
                              fweight=w_user if (weights is not None and weight_type == "fweight") else None)
            num_sep_relu = int(rows.size - m.sum())
            if verbose > 0:
                print(f"(ReLU method dropped {sep_idx.size} separated observations; "
                      f"{num_sep_relu} in all with new singletons)")
            _subset(m)
            y, X, Z, w_user, off_orig = _load()
            solver = FESolver(codes, nlev)
            coll_part = _collinearity(solver)
        stats["t_relu"] += time.perf_counter() - t0

    N = int(y.size)
    Kx, Kd, L0 = len(exog), len(endog), len(inst)
    K = Kx + Kd
    has_offset = off_all is not None
    sumw = float(w_user.sum())
    offset_mean = float(w_user @ off_orig) / sumw if has_offset else 0.0
    off = off_orig - offset_mean if has_offset else off_orig
    F_count = solver.F
    is_fw = weights is not None and weight_type == "fweight"
    N_eff = float(sumw) if is_fw else float(N)

    # ---- IRLS-IV --------------------------------------------------------------------------------
    start_inner_tol = 1e-4
    inner_tol = max(start_inner_tol, tol)
    alt_tol = start_inner_tol
    mean_y = float(w_user @ y) / sumw
    mu = np.full(N, mean_y) if guess == "mean" else 0.5 * (y + mean_y)
    mu = np.maximum(mu, np.maximum(0.05 * y, 1e-3))
    eta = np.log(mu) + off
    pos = y > 0
    zero_idx = np.flatnonzero(~pos)
    with np.errstate(divide="ignore", invalid="ignore"):
        ylogy = np.where(pos, y * np.log(np.where(pos, y, 1.0)), 0.0)
    S_w_ylogy = float(w_user @ ylogy)
    S_wy = float(w_user @ y)
    wy = w_user * y
    sep_mask = np.zeros(N, dtype=bool) if do_sep_mu else None
    num_sep_mu = 0

    b = np.zeros(K)
    converged = False
    ok = 0
    deviance = np.nan
    eps = np.nan
    iter_step_halving = False
    num_step_halving = 0
    old_eta = eta.copy()
    data_dm = None
    z_old = None
    ws_XZ = coll_part if flags["warmstart"] else None
    irls_w = None
    resid = None
    X_dm = Z_dm = z_dm = None
    it = 0
    if verbose > 0:
        print(f"IRLS-IV iterations (N = {N}, K = {K}, L = {Kx + L0})")
    for it in range(1, int(maxit) + 1):
        if not (np.all(np.isfinite(mu)) and np.all(np.isfinite(eta))):
            raise IVPPMLError(f"ivppmlhdfe: mu has infinite or missing values on iteration {it} "
                              "(rc 9003, runaway divergence)")
        if it > 10 and K > 0:
            mb, me = float(np.max(np.abs(b))), float(np.max(np.abs(eta)))
            if mb > 1e6 or me > 30:
                raise IVPPMLError(f"ivppmlhdfe: IRLS diverged at iteration {it} (max|b| = {mb:.2e}, "
                                  f"max|eta| = {me:.2e}; rc 9003)")
        z = eta - off - 1.0 + y / mu
        if do_sep_mu and zero_idx.size:
            z[zero_idx] = eta[zero_idx] - off[zero_idx] - 1.0
        irls_w = np.maximum(w_user * mu, 1e-20)

        tp = time.perf_counter()
        solver.update_weights(irls_w)
        if flags["warmstart"] and data_dm is not None:
            inp = data_dm
            inp[:, 0] += z - z_old
        elif ws_XZ is not None:
            inp = _fstack(z, ws_XZ)
            ws_XZ = None
        else:
            inp = _fstack(z, X, Z)
        data_dm, its, cv = solver.partial_out(inp, inner_tol)
        if flags["warmstart"]:
            z_old = z
        stats["t_partial"] += time.perf_counter() - tp
        stats["inner_iters"] += int(np.sum(its))
        stats["inner_allconv"] &= bool(np.all(cv))
        stats["tol_trace"].append(inner_tol)
        stats["it_trace"].append(int(np.sum(its)))

        z_dm = data_dm[:, 0]
        X_dm = data_dm[:, 1:1 + K]
        Z_dm = _fstack(data_dm[:, 1:1 + Kx], data_dm[:, 1 + K:])
        b_old = b
        b, resid, _ = _tsls(z_dm, X_dm, Z_dm, irls_w)
        beta_change = float(np.max(np.abs(b - b_old))) if K else 0.0

        if not iter_step_halving:
            old_eta = eta
        eta = z - resid + off

        if do_sep_mu and it > 1:
            log_septol = np.log(1e-6)
            adj = log_septol + min(float(eta[pos].min()) + 5.0, 0.0) if pos.any() else log_septol
            sep_mask |= (eta <= adj) & ~pos
            n_new = int(sep_mask.sum()) - num_sep_mu
            if n_new > 0:
                num_sep_mu = int(sep_mask.sum())
                if verbose > 0:
                    print(f"(mu-separation: {n_new} obs detected at iteration {it})")

        mu = np.exp(eta)
        if do_sep_mu and num_sep_mu > 0:
            mu[sep_mask] = 0.0
        np.maximum(mu, EPS100, out=mu)

        old_deviance = deviance
        dev = float(w_user @ mu) - S_wy + S_w_ylogy - float(wy @ eta)
        if 2.0 * dev / N < np.finfo(float).eps:
            dev = 0.0
        if abs(dev) < np.finfo(float).eps:
            dev = 0.0
        deviance = max(2.0 * dev, 0.0)

        if it > 1:
            delta_dev = old_deviance - deviance
            if deviance < 0.1 * delta_dev:
                delta_dev = deviance
            denom = max(min(deviance, old_deviance), 0.1)
            eps = abs(delta_dev) / denom
            if eps < tol:
                if inner_tol <= 1.1 * target_inner_tol or F_count == 1:
                    ok += 1
                    converged = True
                else:
                    # v0.9.4 leaves the step-halving flag set here, which loops to maxit when
                    # this happens right after a step-halving; clear it (documented deviation)
                    iter_step_halving = False
            elif delta_dev < 0 and num_step_halving < 2:
                eta = 0.9 * old_eta + 0.1 * eta
                if num_step_halving > 0:
                    eta = np.maximum(eta, -10.0)
                mu = np.maximum(np.exp(eta), EPS100)
                iter_step_halving = True
                ok = 0
                num_step_halving += 1
            else:
                iter_step_halving = False
                num_step_halving = 0
                ok = 0

        if verbose > 0:
            msg = f"Iter {it:3d}:  dev = {deviance:<11.5e}"
            if it > 1:
                msg += f"  eps = {eps:<9.4e}"
            msg += f"  db = {beta_change:<9.4e}  tol = {inner_tol:5.0e}"
            if iter_step_halving:
                msg += "  H"
            if ok:
                msg += "  O"
            print(msg)

        if iter_step_halving:
            deviance = old_deviance
            continue
        if converged:
            break

        if flags["tolpolicy"] == "fereg":
            k_next = it + 1
            fr = min(k_next - 1, 7) / 7.0
            staged = target_inner_tol if k_next >= 8 else \
                float(np.exp((1 - fr) * np.log(start_inner_tol) + fr * np.log(target_inner_tol)))
            if it > 1:
                prog = target_inner_tol if eps <= tol else \
                    min(max(0.1 * eps, target_inner_tol), start_inner_tol)
            else:
                prog = start_inner_tol
            inner_tol = max(min(staged, prog, inner_tol), target_inner_tol)
        elif it > 1 and eps < inner_tol:
            inner_tol = max(min(0.1 * inner_tol, alt_tol), target_inner_tol)
            alt_tol = 10.0 ** (-np.ceil(np.log10(1.0 / max(0.1 * eps, np.finfo(float).eps))))

    if not converged:
        raise IVPPMLConvergenceError(f"ivppmlhdfe: IRLS failed to converge in {maxit} iterations "
                                     f"(eps = {eps:.4e}; rc 430)")
    if verbose > 0:
        print(f"Converged in {it} iterations (tol = {tol:.4e})")

    # ---- final columns for the VCE (v0.9.4) or refreshed at the converged weights (accel R) ---
    if flags["refresh"]:
        tr = time.perf_counter()
        omega = np.maximum(w_user * mu, 1e-20)
        z_f = eta - off - 1.0 + y / mu
        if do_sep_mu and zero_idx.size:
            z_f[zero_idx] = eta[zero_idx] - off[zero_idx] - 1.0
        solver.update_weights(omega)
        if flags["warmstart"]:
            inp = data_dm
            inp[:, 0] += z_f - z_old
        else:
            inp = _fstack(z_f, X, Z)
        data_dm, its, cv = solver.partial_out(inp, inner_tol)
        stats["inner_iters"] += int(np.sum(its))
        stats["inner_allconv"] &= bool(np.all(cv))
        z_dm = data_dm[:, 0]
        X_dm = data_dm[:, 1:1 + K]
        Z_dm = _fstack(data_dm[:, 1:1 + Kx], data_dm[:, 1 + K:])
        resid = z_dm - X_dm @ b
        w_vce = omega
        stats["t_refresh"] = time.perf_counter() - tr
    else:
        w_vce = irls_w
    Zw = Z_dm * w_vce[:, None]
    Pi = _solve_sym(Zw.T @ Z_dm, Zw.T @ X_dm)
    Xhat = Z_dm @ Pi
    w_meat = np.sqrt(w_user) * mu if is_fw else w_vce
    V, N_clust, rank = sandwich_vce(Xhat, X_dm, w_vce, w_meat, resid, vce_type,
                                    clusters=cl_codes, N_eff=N_eff)

    # ---- constant, likelihoods, FE sum ------------------------------------------------------
    eta_no = eta - off_orig
    sw_v = float(w_vce.sum())
    mean_X = (w_vce @ X) / sw_v if K else np.zeros(0)
    b_cons = float(w_vce @ eta_no) / sw_v - float(mean_X @ b)
    ll = float(w_user @ (y * eta - mu - gammaln(y + 1.0)))
    if has_offset and np.sum(np.abs(off_orig)) > 0:
        c0 = np.log(S_wy / float(w_user @ np.exp(off_orig)))
        mu0 = np.exp(c0 + off_orig)
        ll_0 = float(w_user @ (y * (c0 + off_orig) - mu0 - gammaln(y + 1.0)))
    else:
        ll_0 = float(w_user @ (y * np.log(mean_y) - mean_y - gammaln(y + 1.0))) if mean_y > 0 else 0.0
    xb = X @ b + b_cons
    d = eta_no - xb
    d = d - float(w_vce @ d) / sw_v

    # ---- ffirst -----------------------------------------------------------------------------
    ff = None
    df_a = None
    if ffirst and Kd > 0:
        from ._ffirst import absorbed_dof, first_stage_stats
        tf = time.perf_counter()
        df_a = absorbed_dof(codes, nlev, cl_codes if vce_type == "cluster" else None)
        sdof = df_a if df_a > 0 else 1
        psi = w_vce * w_vce / w_user if is_fw else w_vce * w_vce
        ff = first_stage_stats(X_dm[:, :Kx], X_dm[:, Kx:], Z_dm[:, Kx:], w_vce, psi, vce_type,
                               cl_codes, N_eff, sdof, endog)
        stats["t_ffirst"] = time.perf_counter() - tf

    names = exog + endog
    coef = pd.Series(b, index=names, name="coef")
    vcov_df = pd.DataFrame(V, index=names, columns=names)
    se = pd.Series(np.sqrt(np.abs(np.diag(V))), index=names, name="se")
    stats["t_total"] = time.perf_counter() - t_fit0
    stats["outer_iters"] = it
    stats["restarts"] = solver.restarts
    stats["passes"] = solver.passes
    t_part_all = stats["t_partial"] + stats["t_coll"] + stats["t_refresh"]
    stats["partial_share"] = t_part_all / stats["t_total"] if stats["t_total"] > 0 else np.nan

    return IVPPMLResult(
        version=VERSION, depvar=pf.depvar, exog=exog, endog=endog, instruments=inst,
        absorb=absorb, formula=formula, coef=coef, vcov=vcov_df, se=se, intercept=b_cons,
        N=int(round(N_eff)) if is_fw else N, N_rows=N, N_full=N_full, iterations=it,
        converged=converged, deviance=deviance, ll=ll, ll_0=ll_0,
        num_singletons=num_singletons, num_separated=num_sep_relu + num_sep_mu,
        num_sep_relu=num_sep_relu, num_sep_mu=num_sep_mu, rank=rank,
        vce_type=vce_type, cluster_vars=clustvars, N_clust=N_clust,
        weights=weights, weight_type=weight_type if weights is not None else None,
        offset=offset if exposure is None else f"ln({exposure})",
        accel_flags=flags, tol=tol, itol=target_inner_tol, separation=sorted(sep),
        ffirst=ff, df_a=df_a, stats=stats, esample_idx=rows, N_rows_data=N_rows_data, _index=index,
        _eta=eta, _xb=xb, _d=d,
    )
