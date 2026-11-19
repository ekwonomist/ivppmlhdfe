"""First-stage statistics (ffirst), written from ``stata/SPEC_ffirst.md`` (v1.0.0 build).

At convergence, for the linear IV problem at the final IRLS weights omega = w*mu with the fixed
effects partialled out under omega: the statistics are those that
``ivreghdfe z X1 (X2 = Z) [aw = omega], absorb(<same FEs>) <same vce> ffirst`` reports.

* per endogenous regressor: first-stage F (Wald on the excluded instruments, same VCE type as the
  model, ivreg2's small-sample rule), partial R2, Sanderson-Windmeijer (2016) F;
* system: Kleibergen-Paap (2006) rk LM (underidentification) and rk Wald F, Cragg-Donald Wald F.
Kleibergen-Paap is reimplemented from the paper (no GPL code).
"""
from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

__all__ = ["first_stage_stats", "absorbed_dof"]


# ----------------------------------------------------------------------------------------------
# absorbed degrees of freedom (reghdfe's estimate_dof rule, = ivppmlhdfe's e(df_a))
# ----------------------------------------------------------------------------------------------
def _nested(c, k):
    """Is FE c (codes 0..G-1) nested within cluster k (each level of c has one value of k)?"""
    key = c.astype(np.int64) * (int(k.max()) + 1) + k.astype(np.int64)
    return pd.unique(key).size == pd.unique(c).size


def _mobility_groups(a, b):
    """Connected components of the bipartite graph levels(a) - levels(b), edges = observations."""
    ga, gb = int(a.max()) + 1, int(b.max()) + 1
    g = coo_matrix((np.ones(a.size), (a, b + ga)), shape=(ga + gb, ga + gb))
    ncomp, _ = connected_components(g, directed=False)
    # levels that do not occur are isolated nodes; codes are compact, so every node occurs
    return int(ncomp)


def absorbed_dof(codes, nlev, clusters=None):
    F = len(codes)
    M = [0] * F
    nested = [False] * F
    for g in range(F):
        for k in (clusters or []):
            if _nested(codes[g], k):
                M[g] = nlev[g]
                nested[g] = True
                break
    idx = [g for g in range(F) if not nested[g]]
    if len(idx) > 1:
        for g in idx[1:]:
            M[g] = 1
    if any(nested) and len(idx) >= 1:
        for g in idx:
            M[g] = 1
    for a, b in itertools.combinations(idx, 2):
        M[b] = max(M[b], _mobility_groups(codes[a], codes[b]))
    return int(sum(nlev[g] - M[g] for g in range(F)))


# ----------------------------------------------------------------------------------------------
def _codes(cols):
    key = cols[0].astype(np.int64)
    for c in cols[1:]:
        key = key * (int(c.max()) + 1) + c
    codes, uniq = pd.factorize(key, sort=True)
    return codes.astype(np.int64), len(uniq)


def _sqrtm_psd(A):
    ev, Q = np.linalg.eigh(0.5 * (A + A.T))
    return (Q * np.sqrt(np.maximum(ev, 0.0))) @ Q.T


def first_stage_stats(W, Xe, Zx, omega, psi, vce_type, clusters, N, sdof, names):
    """First-stage statistics.

    W : N_obs x Kw included exogenous (partialled); Xe : N_obs x Ke endogenous; Zx : N_obs x Lz
    excluded instruments; omega : cross-product weights; psi : robust-meat weights; clusters: list
    of integer code arrays; N : e(N); sdof : absorbed dof (1 if 0).
    """
    Ke, Kw, Lz = Xe.shape[1], W.shape[1], Zx.shape[1]
    L = Lz + Kw
    R = np.column_stack([Zx, W])
    is_cl = vce_type == "cluster" and clusters
    if is_cl:
        nclu = len(clusters)
        one_way = [_codes([c]) for c in clusters]
        G = min(g for _, g in one_way)
        subsets = []
        for size in range(1, nclu + 1):
            for sub in itertools.combinations(range(nclu), size):
                cg = one_way[sub[0]] if size == 1 else _codes([clusters[j] for j in sub])
                subsets.append((1.0 if size % 2 == 1 else -1.0, cg[0], cg[1]))

    def meat(S):
        if not is_cl:
            return (S * psi[:, None]).T @ S
        Sw = S * omega[:, None]
        M = np.zeros((S.shape[1], S.shape[1]))
        for sign, cg, g in subsets:
            U = np.empty((g, S.shape[1]))
            for j in range(S.shape[1]):
                U[:, j] = np.bincount(cg, weights=Sw[:, j], minlength=g)
            M += sign * (U.T @ U)
        return M

    def to_F(W0, q):
        if is_cl:
            return W0 / q * (N - L - sdof) / (N - 1) * (G - 1) / G, G - 1
        return W0 / q * (N - L - sdof) / N, N - L - sdof

    def wald(Rinv, S, b):
        V0 = Rinv @ meat(S) @ Rinv
        Vz = V0[:Lz, :Lz]
        bz = b[:Lz]
        return float(bz @ np.linalg.solve(Vz, bz))

    Rw = R * omega[:, None]
    RWR = Rw.T @ R
    RWRi = np.linalg.inv(0.5 * (RWR + RWR.T))
    Pi_all = RWRi @ (Rw.T @ Xe)                      # L x Ke
    E_all = Xe - R @ Pi_all
    out = {}
    if Kw:
        Ww = W * omega[:, None]
        WWi = np.linalg.inv(Ww.T @ W)
    for k in range(Ke):
        e = E_all[:, k]
        W0 = wald(RWRi, R * e[:, None], Pi_all[:, k])
        F, df2 = to_F(W0, Lz)
        p = float(stats.f.sf(F, Lz, df2))
        x = Xe[:, k]
        er = x - W @ (WWi @ (Ww.T @ x)) if Kw else x
        pr2 = 1.0 - float(e @ (omega * e)) / float(er @ (omega * er))
        if Ke == 1:
            swf, swdf1, swdf2, swp = F, float(Lz), float(df2), p
        else:
            others = [j for j in range(Ke) if j != k]
            Xhat = R @ Pi_all
            Sfit = np.column_stack([Xhat[:, others], W])
            Sact = np.column_stack([Xe[:, others], W])
            Sfw = Sfit * omega[:, None]
            b1 = np.linalg.solve(Sfw.T @ Sfit, Sfw.T @ x)
            e1 = x - Sact @ b1
            b2 = RWRi @ (Rw.T @ e1)
            e2 = e1 - R @ b2
            q = Lz - Ke + 1
            W0s = wald(RWRi, R * e2[:, None], b2)
            swf, swdf2 = to_F(W0s, q)
            swdf1 = float(q)
            swp = float(stats.f.sf(swf, q, swdf2))
        out[names[k]] = {"F": float(F), "df1": float(Lz), "df2": float(df2), "p": p,
                         "partialR2": pr2, "SWF": float(swf), "SWdf1": float(swdf1),
                         "SWdf2": float(swdf2), "SWp": float(swp)}

    # ---- system statistics ---------------------------------------------------------------------
    if Kw:
        Yb = Xe - W @ (WWi @ (Ww.T @ Xe))
        Zb = Zx - W @ (WWi @ (Ww.T @ Zx))
    else:
        Yb, Zb = Xe, Zx
    Zbw = Zb * omega[:, None]
    A = Zbw.T @ Zb
    B = (Yb * omega[:, None]).T @ Yb
    C = Zbw.T @ Yb
    Pi = np.linalg.solve(A, C)                       # Lz x Ke
    LA = np.linalg.cholesky(0.5 * (A + A.T))
    LB = np.linalg.cholesky(0.5 * (B + B.T))
    Theta = LA.T @ Pi @ np.linalg.inv(LB.T)
    U, s, Vt = np.linalg.svd(Theta, full_matrices=True)
    V = Vt.T
    kk = Ke - 1
    if kk == 0:
        Aperp = np.eye(Lz)
        Bperp = np.eye(1)
    else:
        u12, u22 = U[:kk, kk:], U[kk:, kk:]
        Aperp = U[:, kk:] @ np.linalg.inv(u22) @ _sqrtm_psd(u22 @ u22.T)
        v22 = V[kk:, kk:]
        Bperp = _sqrtm_psd(v22 @ v22.T) @ np.linalg.inv(v22.T) @ V[:, kk:].T
    T = np.kron(Bperp, Aperp.T)
    lam = T @ Theta.ravel(order="F")
    J = np.kron(np.linalg.inv(LB), np.linalg.inv(LA))
    TJ = T @ J

    def rk(v):
        S = (v[:, :, None] * Zb[:, None, :]).reshape(v.shape[0], Ke * Lz)
        mid = TJ @ meat(S) @ TJ.T
        return float(lam @ np.linalg.solve(0.5 * (mid + mid.T), lam))

    idstat = rk(Yb)
    iddf = Lz - Ke + 1
    rkw = rk(Yb - Zb @ Pi)
    if is_cl:
        widstat = rkw / (N - 1) * (N - L - sdof) * (G - 1) / G / Lz
    else:
        widstat = rkw / N * (N - L - sdof) / Lz
    s2 = s[:Ke] ** 2
    cdf = float(np.min(s2 / (1.0 - s2))) * (N - sdof - L) / Lz
    out.update({"idstat": idstat, "iddf": float(iddf), "idp": float(stats.chi2.sf(idstat, iddf)),
                "cdf": float(cdf), "widstat": float(widstat)})
    return out
