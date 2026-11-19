"""Data preparation: sample, fixed-effect codes, singletons, separation, collinearity.

Conventions of Stata ivppmlhdfe v0.9.4 (``program/stata/ivppmlhdfe.ado``):
* markout of every used variable; y >= 0; weights >= 0 and zero weights dropped;
* fixed-effect stage (reghdfe's ``fixed_effects()`` with the iweight trick of ppmlhdfe): iterated
  over the fixed effects until nothing changes, drop observations in groups with one
  observation and -- with separation "fe" -- in groups whose y is all zero. All of these are
  counted in ``num_singletons`` (reghdfe adds both kinds to ``HDFE.num_singletons``; v0.9.4's
  ``num_sep_fe`` is therefore always 0). ``keep_singletons`` switches the whole stage off, as in
  Stata (``drop_singletons = 0`` also skips the zero-group drop);
* collinearity: exog, endog and instruments partialled out under the user weights at the target
  inner tolerance, then a within-block rank test per block;
* ReLU separation (Correia, Guimaraes and Zylkin 2019, as in ppmlhdfe's ``relu_fix_separation``):
  regressors = exog + endog; dropped observations, plus singletons that the drop creates, are
  counted in ``num_separated``.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import scipy.linalg

__all__ = ["fe_codes", "drop_fe_stage", "compact_codes", "select_not_collinear",
           "relu_separation"]


def _factorize(values):
    codes, uniq = pd.factorize(values, sort=False, use_na_sentinel=True)
    return codes.astype(np.int64), len(uniq)


def fe_codes(df, terms, rows):
    """Integer codes (0..G-1) for each FE term (a list of columns, interacted)."""
    out, nlev = [], []
    for term in terms:
        c, g = _factorize(df[term[0]].to_numpy()[rows])
        for col in term[1:]:
            c2, g2 = _factorize(df[col].to_numpy()[rows])
            c, g = _factorize(c * np.int64(g2) + c2)
        out.append(c)
        nlev.append(g)
    return out, nlev


def compact_codes(codes, keep):
    out, nlev = [], []
    for c in codes:
        cc, g = _factorize(c[keep])
        out.append(cc)
        nlev.append(g)
    return out, nlev


def drop_fe_stage(codes, nlevels, keep, y=None, fweight=None):
    """Iterated drop of singleton groups (and, if ``y`` is given, of groups with all y == 0).

    ``fweight``: if given, a singleton is a group whose fweight sum is 1 (reghdfe's rule when it
    reloads with frequency weights); else a group with one row.
    Returns the updated ``keep`` mask (a copy).
    """
    keep = keep.copy()
    F = len(codes)
    if F == 0:
        return keep
    unchanged = 0
    f = 0
    # cycle over the FEs until F consecutive FEs drop nothing (reghdfe's loop)
    while unchanged < F:
        c, G = codes[f], nlevels[f]
        kf = keep.astype(np.float64)
        if fweight is not None:
            cnt_s = np.bincount(c, weights=kf * fweight, minlength=G)
            bad_g = cnt_s == 1
        else:
            cnt = np.bincount(c, weights=kf, minlength=G)
            bad_g = cnt == 1
        if y is not None:
            cnt = np.bincount(c, weights=kf, minlength=G)
            sy = np.bincount(c, weights=kf * y, minlength=G)
            bad_g |= (cnt > 0) & (sy == 0)
        bad = keep & bad_g[c]
        if bad.any():
            keep &= ~bad
            unchanged = 1 if F > 1 else F
        else:
            unchanged += 1
        f = (f + 1) % F
    return keep


def select_not_collinear(A_partial, A_raw, w, rtol=1e-9):
    """Indices of the columns kept, in order: (a) a column whose FE-partialled weighted norm is
    below rtol * the largest raw weighted norm is absorbed by the fixed effects; (b) among the
    rest, a column is dropped when it lies in the span of the columns kept before it (relative
    residual norm <= rtol), so the earlier variable is kept, as Stata's rmcoll does."""
    K = A_partial.shape[1]
    if K == 0:
        return []
    sw = np.sqrt(np.maximum(w, 0.0))
    raw = np.sqrt(np.einsum("ij,ij,i->j", A_raw, A_raw, sw * sw))
    mx = raw.max()
    if mx == 0:
        return []
    WA = sw[:, None] * A_partial
    part = np.linalg.norm(WA, axis=0)
    nz = [j for j in range(K) if part[j] > rtol * mx]
    if len(nz) <= 1:
        return nz
    # order-preserving rank test: non-pivoted QR, R_jj = residual norm of column j on the
    # columns before it
    R = scipy.linalg.qr(WA[:, nz], mode="r")[0]
    keep = []
    for i, j in enumerate(nz):
        if not keep:
            keep.append(j)
            continue
        # residual of column j on the kept columns (re-solve when a column was dropped before)
        cols = [nz.index(k) for k in keep]
        if len(cols) == i:
            res = abs(R[i, i])
        else:
            Q = WA[:, keep]
            coef = np.linalg.lstsq(Q, WA[:, j], rcond=None)[0]
            res = np.linalg.norm(WA[:, j] - Q @ coef)
        if res > rtol * part[j]:
            keep.append(j)
    return keep


def relu_separation(solver, y, X, inner_tol_target, tol=1e-4, maxiter=100, verbose=0):
    """ReLU separation test. Returns the indices (into y) of separated observations.

    Method of weighting for least squares with equality constraints: weight 1 on y == 0, weight
    M = ceil(eps^-1/2) on y > 0; u = 1(y == 0); iterate u <- max(xbd, 0) on y == 0 until no
    negative fitted value (separated: fitted > 0) or no negative residual. FE inner tolerance
    max(tol^2, 1e-13). ``solver`` is an FESolver on the current sample (its weights are changed).
    """
    is_b = y == 0
    nb = int(is_b.sum())
    if nb == 0:
        return np.zeros(0, dtype=np.int64)
    regtol = max(tol * tol, 1e-13)
    Mw = float(np.ceil(np.finfo(float).eps ** -0.5))
    w = np.where(is_b, 1.0, Mw)
    solver.update_weights(w)
    Xt, _, _ = solver.partial_out(X, regtol)
    sw = np.sqrt(w)
    A = sw[:, None] * Xt
    norm_unc = np.abs(X[is_b]).sum(axis=0).max() if X.shape[1] else 0.0
    norm_con = np.abs(X[~is_b]).sum(axis=0).max() if X.shape[1] else 0.0
    ratio = max(norm_unc / norm_con, 1.0) if norm_con > 0 else 1.0
    u = is_b.astype(np.float64)
    ut = None
    u_last = None
    sep = None
    for it in range(1, maxiter + 1):
        ut_in = u if it == 1 else ut + (u - u_last)
        u_last = u.copy()
        ut, _, _ = solver.partial_out(ut_in, regtol)
        if X.shape[1]:
            bcoef = np.linalg.lstsq(A, sw * ut, rcond=None)[0]
            resid = ut - Xt @ bcoef
        else:
            resid = ut.copy()
        yy = float(u @ u)
        epsilon = max(yy / Mw * ratio ** 2, 1e-8)
        delta = epsilon + tol
        xbd = u - resid
        xbd[~is_b] = 0.0
        xb = xbd[is_b]
        xb[(xb >= -0.1 * delta) & (xb <= delta)] = 0.0
        if verbose > 1:
            print(f"  relu iter {it}: u'u={yy:.6g} min(xbd)={xb.min():.3g} candidates={(xb > delta).sum()}")
        if np.all(xb >= 0):
            sep = xb > 0
            break
        resid[np.abs(resid) < 1e-8] = 0.0
        rb = resid[is_b]
        if rb.min() >= 0:
            xb[rb > delta] = 0.0
            sep = xb > 0
            break
        if it == maxiter:
            warnings.warn("ivppmlhdfe: ReLU separation check reached its maximum number of "
                          "iterations; no observation dropped", RuntimeWarning, stacklevel=3)
            return np.zeros(0, dtype=np.int64)
        u[is_b] = np.maximum(xb, 0.0)
    idx_b = np.flatnonzero(is_b)
    return idx_b[sep]
