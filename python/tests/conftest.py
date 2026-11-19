import json
import os
import sys
from pathlib import Path

import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "src"
if not os.environ.get("IVPPML_INSTALLED"):        # set it to test the pip-installed package
    sys.path.insert(0, str(SRC))

# testdata: IVPPML_TESTDATA, else the repository's crosscheck/ folder (crosscheck/testdata, crosscheck/results)
_CROSS = HERE.parent.parent / "crosscheck"
TESTDATA = Path(os.environ.get("IVPPML_TESTDATA", _CROSS / "testdata"))
RESULTS = Path(os.environ.get("IVPPML_RESULTS", _CROSS / "results"))


@pytest.fixture(scope="session")
def cases():
    f = TESTDATA / "cases.json"
    if not f.exists():
        pytest.skip(f"no cases.json at {f}")
    return json.loads(f.read_text(encoding="utf-8"))["cases"]


@pytest.fixture(scope="session")
def datasets():
    cache = {}

    def get(name):
        if name not in cache:
            cache[name] = pd.read_csv(TESTDATA / name)
        return cache[name]
    return get


def case_args(c):
    fe = " + ".join("^".join(t) for t in c["absorb"])
    rhs = " + ".join(c["exog"]) if c["exog"] else "1"
    f = f"{c['depvar']} ~ {rhs}"
    if fe:
        f += f" | {fe}"
    if c["endog"]:
        f += f" | {' + '.join(c['endog'])} ~ {' + '.join(c['instruments'])}"
    v = c["vce"]
    vcov = "robust" if v["type"] == "robust" else {"CL": list(v["vars"])}
    kw = dict(vcov=vcov, weights=c.get("weights"), offset=c.get("offset"))
    if c.get("weight_type"):
        kw["weight_type"] = c["weight_type"]
    return f, kw
