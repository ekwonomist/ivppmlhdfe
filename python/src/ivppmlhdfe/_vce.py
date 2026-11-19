"""Sandwich VCE for IV-PPML, with the conventions of Stata ivppmlhdfe v0.9.4.

bread = (Xhat' W X~)^-1, scores_i = Xhat_i * w_meat_i * resid_i.
* robust: V = N/(N-1) bread meat bread, N = sum of fweights with fweight, else the row count;
  meat weights: irls_w = w*mu, except fweight where w_meat = sqrt(w)*mu (the meat then scales
  with w, as if rows were replicated).
* one-way cluster: V = G/(G-1) bread meat bread, scores with irls_w (all weight types);
  G = 1 falls back to robust.
* multi-way cluster (Cameron-Gelbach-Miller): sum over the non-empty subsets of the cluster
  variables of (-1)^(|s|+1) G_s/(G_s-1) bread meat_s bread (subsets with G_s = 1 skipped), then
  symmetrised, and if not PSD, V - min_eig * I (reghdfe's fix). N_clust = min of the one-way G.
"""
from __future__ import annotations

import itertools

import numpy as np
import pandas as pd

__all__ = ["sandwich_vce", "cluster_sums"]


def _inv_sym(A):
    A = 0.5 * (A + A.T)
    try:
        return np.linalg.inv(A)
    except np.linalg.LinAlgError:
        return np.linalg.pinv(A)


def cluster_sums(scores, codes, G):
    """Per-cluster column sums of scores (N x K) for integer codes 0..G-1."""
    K = scores.shape[1]
    out = np.empty((G, K))
    for k in range(K):
        out[:, k] = np.bincount(codes, weights=scores[:, k], minlength=G)
    return out


def _codes_of(cols):
    """Compact codes for the interaction of one or more integer code arrays."""
    if len(cols) == 1:
        codes, uniq = pd.factorize(cols[0], sort=True)
        return codes.astype(np.int64), len(uniq)
    key = cols[0].astype(np.int64)
    for c in cols[1:]:
        key = key * (int(c.max()) + 1) + c
    codes, uniq = pd.factorize(key, sort=True)
    return codes.astype(np.int64), len(uniq)


def sandwich_vce(Xhat, X_dm, w_bread, w_meat, resid, vce_type, clusters=None, N_eff=None):
    """Return (V, N_clust list, actual_rank).

    clusters: list of integer code arrays (0..G-1), one per cluster variable.
    """
    N, K = X_dm.shape
    if K == 0:
        return np.zeros((0, 0)), [], 0
    bread_inv = Xhat.T @ (w_bread[:, None] * X_dm)
    bread = _inv_sym(bread_inv)
    rank = int(np.linalg.matrix_rank(0.5 * (bread_inv + bread_inv.T)))
    if N_eff is None:
        N_eff = N
    if vce_type != "cluster":
        sc = Xhat * (w_meat * resid)[:, None]
        meat = sc.T @ sc
        V = (N_eff / (N_eff - 1.0)) * bread @ meat @ bread
        return 0.5 * (V + V.T), [], rank
    sc = Xhat * (w_bread * resid)[:, None]
    nclu = len(clusters)
    if nclu == 1:
        codes, G = _codes_of([clusters[0]])
        if G <= 1:
            sc_r = Xhat * (w_meat * resid)[:, None]
            V = (N_eff / (N_eff - 1.0)) * bread @ (sc_r.T @ sc_r) @ bread
            return 0.5 * (V + V.T), [G], rank
        S = cluster_sums(sc, codes, G)
        V = (G / (G - 1.0)) * bread @ (S.T @ S) @ bread
        return 0.5 * (V + V.T), [G], rank
    V = np.zeros((K, K))
    Gs = [0] * nclu
    for size in range(1, nclu + 1):
        for subset in itertools.combinations(range(nclu), size):
            codes, G = _codes_of([clusters[j] for j in subset])
            if size == 1:
                Gs[subset[0]] = G
            if G > 1:
                S = cluster_sums(sc, codes, G)
                sign = 1.0 if size % 2 == 1 else -1.0
                V += sign * (G / (G - 1.0)) * bread @ (S.T @ S) @ bread
    V = 0.5 * (V + V.T)
    ev = np.linalg.eigvalsh(V)
    if ev.min() < 0:
        V = V - ev.min() * np.eye(K)
    return V, Gs, rank
