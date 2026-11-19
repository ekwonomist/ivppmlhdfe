"""ivppmlhdfe: IV-PPML with high-dimensional fixed effects.

    import ivppmlhdfe
    r = ivppmlhdfe.fit(df, "y ~ x2 | exp^year + imp^year | x1 ~ z", vcov={"CL": ["pair"]})
    print(r.summary())

Authors: Ohyun Kwon, Mario Larch, Jangsu Yoon, Yoto V. Yotov. License: MIT.
"""
from ._fe import FESolver
from ._fit import VERSION as __version__
from ._fit import IVPPMLConvergenceError, IVPPMLError, fit
from ._formula import parse_formula
from ._result import IVPPMLResult

__all__ = ["fit", "IVPPMLResult", "IVPPMLError", "IVPPMLConvergenceError", "FESolver",
           "parse_formula", "__version__"]
