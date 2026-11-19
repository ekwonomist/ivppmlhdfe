"""Fixed-effect engine: exact within transform, dense least-squares checks, warm-start identity."""
import numpy as np
import pytest

from ivppmlhdfe import FESolver


def _dense_resid(V, codes, nlev, w):
    """Weighted LS residuals of V on all FE dummies (dense, small N)."""
    D = np.concatenate([np.eye(g)[c] for c, g in zip(codes, nlev)], axis=1)
    sw = np.sqrt(w)
    coef, *_ = np.linalg.lstsq(sw[:, None] * D, sw[:, None] * V, rcond=None)
    return V - D @ coef


def _codes(rng, N, G):
    c = rng.integers(0, G, N)
    c[:G] = np.arange(G)          # every level present
    return c


def test_one_fe_is_exact_within_transform():
    rng = np.random.default_rng(1)
    N, G = 500, 37
    c = _codes(rng, N, G)
    w = rng.uniform(0.1, 3.0, N)
    V = rng.standard_normal((N, 3))
    S = FESolver([c], [G]).update_weights(w)
    E, its, conv = S.partial_out(V, tol=1e-8)
    means = np.stack([np.bincount(c, weights=w * V[:, j], minlength=G) /
                      np.bincount(c, weights=w, minlength=G) for j in range(3)], axis=1)
    assert np.max(np.abs(E - (V - means[c]))) < 1e-13
    assert conv.all() and its.sum() == 0


@pytest.mark.parametrize("F", [2, 3])
def test_multi_fe_matches_dense_lstsq(F):
    rng = np.random.default_rng(10 + F)
    N = 400
    nlev = [25, 12, 8][:F]
    codes = [_codes(rng, N, g) for g in nlev]
    w = rng.uniform(0.05, 5.0, N)
    V = rng.standard_normal((N, 4))
    S = FESolver(codes, nlev).update_weights(w)
    E, its, conv = S.partial_out(V, tol=1e-13)
    Ed = _dense_resid(V, codes, nlev, w)
    assert conv.all()
    assert np.max(np.abs(E - Ed)) < 1e-9 * np.max(np.abs(V))


def test_disconnected_graph():
    # two separate blocks of (id, year): the Schur system is singular but consistent
    rng = np.random.default_rng(3)
    a = np.repeat(np.arange(10), 6)              # 10 ids x 6 obs
    b = np.where(a < 5, rng.integers(0, 4, a.size), 4 + rng.integers(0, 3, a.size))
    b[:4] = np.arange(4)
    b[30:33] = 4 + np.arange(3)
    w = rng.uniform(0.5, 2, a.size)
    V = rng.standard_normal((a.size, 2))
    S = FESolver([a, b], [10, 7]).update_weights(w)
    E, _, conv = S.partial_out(V, tol=1e-13)
    assert conv.all()
    assert np.max(np.abs(E - _dense_resid(V, [a, b], [10, 7], w))) < 1e-9


def test_warm_start_identity():
    # M_w2 (M_w1 v) = M_w2 v: partialling a previous residual under new weights is exact
    rng = np.random.default_rng(4)
    N = 600
    nlev = [30, 20, 10]
    codes = [_codes(rng, N, g) for g in nlev]
    V = rng.standard_normal((N, 2))
    S = FESolver(codes, nlev)
    E1, _, _ = S.update_weights(rng.uniform(0.2, 2, N)).partial_out(V, tol=1e-13)
    w2 = rng.uniform(0.2, 2, N)
    S.update_weights(w2)
    Ea, _, _ = S.partial_out(V, tol=1e-13)
    Eb, _, _ = S.partial_out(E1, tol=1e-13)
    assert np.max(np.abs(Ea - Eb)) < 1e-10


def test_absorbed_column_terminates():
    # a column that is exactly a fixed effect: residual ~ 0, the solver must stop quickly
    rng = np.random.default_rng(5)
    N = 300
    c1, c2 = _codes(rng, N, 20), _codes(rng, N, 9)
    v = rng.standard_normal(9)[c2]
    S = FESolver([c1, c2], [20, 9]).update_weights(rng.uniform(0.5, 2, N))
    E, its, conv = S.partial_out(np.column_stack([v, rng.standard_normal(N)]), tol=1e-9,
                                 maxiter=2000)
    assert np.max(np.abs(E[:, 0])) < 1e-6
    assert its.max() < 2000
