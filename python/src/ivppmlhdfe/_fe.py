"""Fixed-effect engine: exact elimination of the largest FE + Jacobi-PCG on the rest.

Idea (S) of SPEC section 2, after fereg (Schmieder 2026, NBER WP 36017, pp.22-24) and the
reference implementation ``IVPPMLFast/src/schur.jl``:

* Observations are sorted once by the fixed effect with the most levels (FE 1). Its weighted
  within-group demeaning M1 is exact.
* The remaining FEs (design D2, ``m`` levels in total) are solved from the Schur complement
  S = D2' W M1 D2 by conjugate gradients with a Jacobi preconditioner (the exact diagonal of
  S). S is never formed: S p = D2' W M1 D2 p, one pass over the data for all columns at once.
* Stopping rule per column: ||A'e||_2 <= tol * sqrt(F) * ||sqrt(W) e||_2, with A the
  column-normalised design (so the FE-1 block of A'e is zero by construction and the rest is
  (D2' W e) / sqrt(d2)). A column that passes on the recursive CG residual must pass again on a
  freshly computed residual (fereg's certificate), else CG restarts from the current iterate.
* With one FE the result is the exact weighted within transform. Disconnected graphs make S
  singular but consistent; CG on the projection is fine.

Layout: internally the columns are stored as rows of a C-ordered (nc x N) array in the sorted
order, so every per-observation operation runs over contiguous memory; the gather is ``np.take``
and the scatter ``np.bincount``. No Python loops over observations.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = ["FESolver"]


class FESolver:
    """Partial fixed effects out of the columns of a data matrix under observation weights.

    Parameters
    ----------
    codes : list of int arrays
        One array per fixed effect, each of length N with values 0..G_f-1 and every level present.
    nlevels : list of int
        Number of levels of each fixed effect.

    The structure (sort order, level indices) is built once; ``update_weights`` refills only the
    weight-dependent pieces, so one solver serves every IRLS iteration.
    """

    def __init__(self, codes, nlevels):
        if len(codes) == 0:
            raise ValueError("FESolver needs at least one fixed effect")
        self.N = N = int(len(codes[0]))
        self.F = F = len(codes)
        nlevels = [int(g) for g in nlevels]
        big = int(np.argmax(nlevels))           # ties: the first
        self.big = big
        g = np.asarray(codes[big], dtype=np.intp)
        self.G = G = nlevels[big]
        perm = np.argsort(g, kind="stable")
        self.isid = bool(np.array_equal(perm, np.arange(N)))
        self.perm = perm
        gs = g[perm]
        counts = np.bincount(gs, minlength=G)
        if np.any(counts == 0):
            raise ValueError("FESolver: every level of the eliminated FE must be present")
        self.counts = counts
        self.starts = np.concatenate(([0], np.cumsum(counts)[:-1])).astype(np.intp)
        self.gidx = gs
        # balanced eliminated FE (all groups the same size): demeaning by a reshape
        self.gsize = int(counts[0]) if np.all(counts == counts[0]) else 0
        self._bufs = {}
        rest = [f for f in range(F) if f != big]
        self.R = R = len(rest)
        self.sqrtF = float(np.sqrt(F))
        self.m = 0
        self.levs = []
        if R:
            offs = np.cumsum([0] + [nlevels[f] for f in rest])
            self.levs = [np.asarray(codes[f], dtype=np.intp)[perm] + int(offs[r])
                         for r, f in enumerate(rest)]
            m = int(offs[-1])
            self.m = m
            # unique (eliminated group, remaining level) pairs, for the exact Jacobi diagonal
            key = np.concatenate([gs.astype(np.int64) * np.int64(m) + lv for lv in self.levs])
            inv, uniq = pd.factorize(key, sort=False)
            self.pair_inv = inv.astype(np.intp)
            uniq = np.asarray(uniq, dtype=np.int64)
            self.pair_g = (uniq // m).astype(np.intp)
            self.pair_l = (uniq % m).astype(np.intp)
            self.lev_all = np.concatenate(self.levs)
        # counters (cumulative, reset by the caller)
        self.passes = 0
        self.restarts = 0
        self.w = None

    # ------------------------------------------------------------------------------------------
    def update_weights(self, w):
        """Refill the weight-dependent pieces (group sums, level sums, Jacobi diagonal of S)."""
        w = np.asarray(w, dtype=np.float64)
        ws = np.ascontiguousarray(w if self.isid else w[self.perm])
        self.w = ws
        n1 = np.add.reduceat(ws, self.starts)
        self.n1inv = np.where(n1 > 0, 1.0 / np.where(n1 > 0, n1, 1.0), 0.0)
        if self.R:
            m = self.m
            wr = np.tile(ws, self.R)
            d2 = np.bincount(self.lev_all, weights=wr, minlength=m)
            C = np.bincount(self.pair_inv, weights=wr)
            dS = d2 - np.bincount(self.pair_l, weights=C * C * self.n1inv[self.pair_g], minlength=m)
            ok = dS > 1e-12 * d2
            self.pinv = np.where(ok, 1.0 / np.where(ok, dS, 1.0), 0.0)
            self.d2inv = np.where(d2 > 0, 1.0 / np.where(d2 > 0, d2, 1.0), 0.0)
        return self

    # ------------------------------------------------------------------------------------------
    # kernels on (nc x N) arrays in the sorted order; U and T are reusable work buffers
    def _buffers(self, nc):
        b = self._bufs.get(nc)
        if b is None:
            b = (np.empty((nc, self.N)), np.empty((nc, self.N)))
            self._bufs = {nc: b}            # keep one size only
        return b

    def _within(self, U, T):
        """In place: U <- M1 U (weighted within-group demeaning for the eliminated FE)."""
        np.multiply(U, self.w, out=T)
        if self.gsize:
            nc, G, s = U.shape[0], self.G, self.gsize
            sb = T.reshape(nc, G, s).sum(axis=2)
            sb *= self.n1inv
            U.reshape(nc, G, s)[...] -= sb[:, :, None]
        else:
            sb = np.add.reduceat(T, self.starts, axis=1)
            sb *= self.n1inv
            np.take(sb, self.gidx, axis=1, out=T, mode="clip")
            U -= T
        return U

    def _gather(self, P, U, T):
        # mode="clip": with out= and the default mode="raise" numpy buffers the output; the
        # indices are valid, so clipping never changes the result
        np.take(P, self.levs[0], axis=1, out=U, mode="clip")
        for lv in self.levs[1:]:
            np.take(P, lv, axis=1, out=T, mode="clip")
            U += T
        return U

    def _scatter(self, wE):
        nc = wE.shape[0]
        Q = np.empty((nc, self.m))
        for c in range(nc):
            q = np.bincount(self.levs[0], weights=wE[c], minlength=self.m)
            for lv in self.levs[1:]:
                q += np.bincount(lv, weights=wE[c], minlength=self.m)
            Q[c] = q
        return Q

    def _resid(self, V, X, U, T):
        """U <- E = M1 (V - D2 X); returns (E, Gr = D2' W E, en2 = ||sqrt(W) E||^2 per column)."""
        if X is None:
            np.copyto(U, V)
        else:
            self._gather(X, U, T)
            np.subtract(V, U, out=U)
        self._within(U, T)
        np.multiply(U, self.w, out=T)
        en2 = np.einsum("ij,ij->i", T, U)
        Gr = self._scatter(T)
        self.passes += 1
        return U, Gr, en2

    def _apply(self, P, U, T):
        """Q = S P = D2' W M1 D2 P, all columns in one pass."""
        self._gather(P, U, T)
        self._within(U, T)
        U *= self.w
        self.passes += 1
        return self._scatter(U)

    def _gs(self, A):
        return np.einsum("ij,ij,j->i", A, A, self.d2inv)

    # ------------------------------------------------------------------------------------------
    def partial_out(self, data, tol=1e-8, maxiter=10000):
        """Return the residuals of every column of ``data`` (N x nc) on the fixed effects.

        Returns (residuals as an N x nc Fortran-ordered array in the original row order, CG
        iterations per column, converged per column). A 1-d input gives a 1-d output.
        """
        data = np.asarray(data, dtype=np.float64)
        squeeze = data.ndim == 1
        if squeeze:
            data = data[:, None]
        nc = data.shape[1]
        iters = np.zeros(nc, dtype=np.int64)
        if self.isid:
            V = np.ascontiguousarray(data.T)
            if self.R == 0 and np.shares_memory(V, data):
                V = V.copy()                    # _within works in place
        else:
            V = np.take(data.T, self.perm, axis=1)
        if self.R == 0 or nc == 0:
            E = self._within(V, np.empty_like(V))
            self.passes += 1
            out = self._unsort(E)
            return (out[:, 0] if squeeze else out), iters, np.ones(nc, dtype=bool)

        U, T = self._buffers(nc)
        thr2 = (float(tol) * self.sqrtF) ** 2
        X = np.zeros((nc, self.m))
        _, B, c0 = self._resid(V, None, U, T)
        # floor for (nearly) fully absorbed columns: below 1e-7 relative to the input, the test
        # becomes absolute, else rounding noise could keep CG running forever
        vv = np.einsum("ij,ij,j->i", V, V, self.w)
        floor = 1e-14 * vv
        Rr = B.copy()
        conv = np.zeros(nc, dtype=bool)
        active = self._gs(Rr) > thr2 * np.maximum(c0, floor)
        rz = np.zeros(nc)
        total = 0
        restarts = 0
        pinv = self.pinv
        while True:
            Z = Rr * pinv
            rz = np.where(active, np.einsum("ij,ij->i", Rr, Z), rz)
            active &= rz > 0
            P = np.where(active[:, None], Z, 0.0)
            while active.any() and total < maxiter:
                Q = self._apply(P, U, T)
                total += 1
                pq = np.einsum("ij,ij->i", P, Q)
                active &= pq > 0
                alpha = np.where(active, rz / np.where(active, pq, 1.0), 0.0)[:, None]
                X += alpha * P
                Rr -= alpha * Q
                iters += active
                xb = np.einsum("ij,ij->i", X, B + Rr)
                e2 = np.maximum(c0 - xb, 0.0)
                done = active & (self._gs(Rr) <= thr2 * np.maximum(e2, floor))
                active &= ~done
                Z = Rr * pinv
                s = np.einsum("ij,ij->i", Rr, Z)
                beta = np.where(active, s / np.where(rz != 0, rz, 1.0), 0.0)
                rz = np.where(active, s, rz)
                P = np.where(active[:, None], Z + beta[:, None] * P, 0.0)
            # certificate: recompute the residual and the gradient from scratch
            E, Gr, en2 = self._resid(V, X, U, T)
            ok = self._gs(Gr) <= thr2 * np.maximum(en2, floor)
            conv |= ok
            retry = ~conv & (total < maxiter) & (restarts < 50)
            if not retry.any():
                break
            Rr[retry] = Gr[retry]
            active = retry.copy()
            restarts += 1
        self.restarts += restarts
        out = self._unsort(E)
        return (out[:, 0] if squeeze else out), iters, conv

    def _unsort(self, E):
        """(nc x N) sorted -> (N x nc) Fortran-ordered in the original order (always a fresh
        array: E may be a work buffer)."""
        if self.isid:
            return np.array(E.T, order="F", copy=True)
        out = np.empty_like(E)
        out[:, self.perm] = E
        return out.T
