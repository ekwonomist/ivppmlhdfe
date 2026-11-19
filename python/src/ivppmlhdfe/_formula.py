"""Formula parsing (fixest style).

    "y ~ x2 + x4 | id + year | x1 + x3 ~ z + z2 + z3"

* part 1: ``depvar ~ exogenous regressors`` (``1`` or nothing for none; ``0``/``-1`` accepted too:
  the constant is always absorbed by the fixed effects or by the implicit intercept);
* optional fixed-effect part: terms joined by ``+``; an interaction is written ``exp^year``;
  ``0`` or ``1`` means no fixed effects;
* optional IV part (the last part, containing ``~``): ``endogenous ~ excluded instruments``.
No IV part = plain PPML. No fixed-effect part = an intercept only (Stata's ``noabsorb``).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

__all__ = ["ParsedFormula", "parse_formula"]

_NAME = re.compile(r"^[A-Za-z_.][A-Za-z0-9_.]*$")


@dataclass
class ParsedFormula:
    depvar: str
    exog: list = field(default_factory=list)
    endog: list = field(default_factory=list)
    instruments: list = field(default_factory=list)
    absorb: list = field(default_factory=list)      # list of lists of column names

    def columns(self):
        cols = [self.depvar] + self.exog + self.endog + self.instruments
        for fe in self.absorb:
            cols += fe
        return list(dict.fromkeys(cols))


def _terms(s, what):
    s = s.strip()
    if s in ("", "1", "0", "-1", "+1", "1+0", "0+1"):
        return []
    out = []
    for t in s.split("+"):
        t = t.strip()
        if t in ("", "1", "0", "-1"):
            continue
        if not _NAME.match(t):
            raise ValueError(f"ivppmlhdfe: cannot parse {what} term {t!r} (use plain column names)")
        out.append(t)
    return out


def _fe_terms(s):
    s = s.strip()
    if s in ("", "0", "1"):
        return []
    out = []
    for t in s.split("+"):
        t = t.strip()
        if t in ("", "0", "1"):
            continue
        parts = [p.strip() for p in re.split(r"\^|#", t)]
        for p in parts:
            if not _NAME.match(p):
                raise ValueError(f"ivppmlhdfe: cannot parse fixed-effect term {t!r}")
        out.append(parts)
    return out


def parse_formula(formula: str) -> ParsedFormula:
    if not isinstance(formula, str):
        raise TypeError("ivppmlhdfe: formula must be a string such as 'y ~ x2 | id + year | x1 ~ z'")
    parts = [p.strip() for p in formula.split("|")]
    if "~" not in parts[0]:
        raise ValueError("ivppmlhdfe: the first part of the formula must be 'depvar ~ regressors'")
    lhs, rhs = parts[0].split("~", 1)
    depvar = lhs.strip()
    if not _NAME.match(depvar):
        raise ValueError(f"ivppmlhdfe: cannot parse the dependent variable {depvar!r}")
    exog = _terms(rhs, "regressor")
    absorb, endog, inst = [], [], []
    rest = parts[1:]
    if rest and "~" in rest[-1]:
        iv = rest.pop()
        e, z = iv.split("~", 1)
        endog = _terms(e, "endogenous")
        inst = _terms(z, "instrument")
        if not endog:
            raise ValueError("ivppmlhdfe: the IV part needs at least one endogenous regressor")
    if len(rest) > 1:
        raise ValueError("ivppmlhdfe: formula has too many '|' parts")
    if rest:
        if "~" in rest[0]:
            raise ValueError("ivppmlhdfe: the IV part must be the last part of the formula")
        absorb = _fe_terms(rest[0])
    for lst, nm in ((exog, "regressor"), (endog, "endogenous regressor"), (inst, "instrument")):
        if len(set(lst)) != len(lst):
            raise ValueError(f"ivppmlhdfe: a {nm} appears twice")
    return ParsedFormula(depvar, exog, endog, inst, absorb)
