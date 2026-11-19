"""Result object of ivppmlhdfe.fit."""
from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy import stats

__all__ = ["IVPPMLResult"]


class IVPPMLResult:
    """Estimates of an IV-PPML model with high-dimensional fixed effects.

    Main fields: ``coef`` (pandas Series, slopes only), ``vcov`` (DataFrame), ``se`` (Series),
    ``intercept`` (the partialled-out constant, reported without an SE as in Stata), ``N``
    (sum of fweights with fweight, else rows), ``N_full``, ``iterations``, ``converged``,
    ``deviance``, ``ll``, ``ll_0``, ``num_singletons``, ``num_separated``, ``N_clust``,
    ``ffirst`` (first-stage statistics or None), ``stats`` (timings and inner iterations).
    """

    def __init__(self, **kw):
        self.__dict__.update(kw)

    # -- accessors ------------------------------------------------------------------------------
    @property
    def coefnames(self):
        return list(self.coef.index)

    def nobs(self):
        return self.N

    def tstat(self):
        return self.coef / self.se

    def pvalues(self):
        z = self.tstat()
        return pd.Series(2 * stats.norm.sf(np.abs(z.to_numpy())), index=self.coef.index)

    def confint(self, level=0.95):
        q = stats.norm.ppf(0.5 + level / 2)
        return pd.DataFrame({"lower": self.coef - q * self.se, "upper": self.coef + q * self.se})

    @property
    def pseudo_r2(self):
        return 1 - self.ll / self.ll_0 if self.ll_0 else float("nan")

    # -- display --------------------------------------------------------------------------------
    def summary(self, digits=6):
        lines = []
        title = "IV-PPML with high-dimensional fixed effects" if self.endog else \
            "PPML with high-dimensional fixed effects"
        lines.append(f"ivppmlhdfe {self.version}: {title}")
        lines.append(f"Dependent variable: {self.depvar}")
        absorb = " ".join("#".join(t) for t in self.absorb) if self.absorb else "_cons (none)"
        lines.append(f"Absorbed FE: {absorb}")
        if self.endog:
            lines.append(f"Endogenous: {' '.join(self.endog)}   Excluded instruments: "
                         f"{' '.join(self.instruments)}")
        vce = "robust" if self.vce_type != "cluster" else \
            "cluster " + " ".join(self.cluster_vars)
        lines.append(f"N = {self.N:,}   iterations = {self.iterations}   "
                     f"deviance = {self.deviance:.6g}   ll = {self.ll:.6g}   "
                     f"pseudo-R2 = {self.pseudo_r2:.4f}")
        extra = f"VCE: {vce}"
        if self.vce_type == "cluster":
            extra += f"   N_clust = {min(self.N_clust):,}" if self.N_clust else ""
        extra += f"   singletons dropped = {self.num_singletons}   separated = {self.num_separated}"
        lines.append(extra)
        acc = "on" if self.accel_flags == {"warmstart": True, "tolpolicy": "fereg", "refresh": True} \
            else ("off" if not any(v not in (False, "ppmlhdfe") for v in self.accel_flags.values())
                  else str(self.accel_flags))
        lines.append(f"accel: {acc}")
        hdr = f"{'':>14} {'Coef.':>13} {'Std. Err.':>13} {'z':>9} {'P>|z|':>8} {'[95% Conf. Interval]':>27}"
        lines.append("-" * len(hdr))
        lines.append(hdr)
        lines.append("-" * len(hdr))
        ci = self.confint()
        p = self.pvalues()
        z = self.tstat()
        for k in self.coef.index:
            lines.append(f"{k:>14} {self.coef[k]:>13.{digits}g} {self.se[k]:>13.{digits}g} "
                         f"{z[k]:>9.2f} {p[k]:>8.3f} {ci.loc[k, 'lower']:>13.{digits}g} "
                         f"{ci.loc[k, 'upper']:>13.{digits}g}")
        lines.append(f"{'_cons':>14} {self.intercept:>13.{digits}g} {'(absorbed)':>13}")
        lines.append("-" * len(hdr))
        if self.ffirst:
            lines.append(_ffirst_text(self.ffirst))
        return "\n".join(lines)

    def __repr__(self):
        return self.summary()

    def __str__(self):
        return self.summary()

    # -- prediction -----------------------------------------------------------------------------
    def predict(self, kind="mu"):
        """In-sample predictions, aligned with the rows of the estimation data (NaN outside the
        estimation sample): ``mu`` = exp(xb + d + offset), ``eta`` = xb + d + offset,
        ``xb`` = X b + _cons, ``d`` = the sum of the fixed effects (weighted mean 0)."""
        if kind not in ("mu", "eta", "xb", "d"):
            raise ValueError("predict: kind must be 'mu', 'eta', 'xb' or 'd'")
        vals = {"mu": np.exp(self._eta), "eta": self._eta, "xb": self._xb, "d": self._d}[kind]
        out = np.full(self.N_rows_data, np.nan)
        out[self.esample_idx] = vals
        if self._index is not None:
            return pd.Series(out, index=self._index, name=kind)
        return out

    @property
    def esample(self):
        m = np.zeros(self.N_rows_data, dtype=bool)
        m[self.esample_idx] = True
        return m


def _ffirst_text(ff):
    lines = ["First-stage statistics (ffirst)"]
    for k, v in ff.items():
        if isinstance(v, dict):
            lines.append(f"  {k}: F({v.get('df1', float('nan')):g}, {v.get('df2', float('nan')):g}) = "
                         f"{v.get('F', float('nan')):.4f} (p = {v.get('p', float('nan')):.4g}), "
                         f"partial R2 = {v.get('partialR2', float('nan')):.4f}, "
                         f"SW F = {v.get('SWF', float('nan')):.4f}")
    if "idstat" in ff:
        lines.append(f"  Kleibergen-Paap rk LM = {ff['idstat']:.4f} (df {ff['iddf']:g}, "
                     f"p = {ff['idp']:.4g});  Cragg-Donald Wald F = {ff['cdf']:.4f};  "
                     f"Kleibergen-Paap rk Wald F = {ff['widstat']:.4f}")
    return "\n".join(lines)


def _isnan(x):
    return x is None or (isinstance(x, float) and math.isnan(x))
