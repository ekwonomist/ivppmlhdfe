"""Compare two ivppmlhdfe result files (the shared JSON schema in SPEC.md, section 7).

    py -3.12 program/sandbox/v010_build/tools/compare.py A.json B.json [--tol 1e-6] [--ffirst]

Relative difference = |a - b| / max(|b|, 1e-3), so a coefficient near zero is compared on an
absolute scale. Counts (N, num_singletons, num_separated) must match exactly.
Exit code 0 if every case present in both files is within --tol, 1 otherwise.
"""
import argparse
import json
import math
import sys


def rel(a, b):
    if a is None or b is None:
        return math.inf
    return abs(a - b) / max(abs(b), 1e-3)


def maxrel(da, db):
    if not da and not db:
        return 0.0, ""
    worst, where = 0.0, ""
    for k in set(da) | set(db):
        r = rel(da.get(k), db.get(k))
        if r > worst:
            worst, where = r, k
    return worst, where


def ffirst_flat(ff):
    out = {}
    if not ff:
        return out
    for k, v in ff.items():
        if isinstance(v, dict):
            for s, x in v.items():
                if isinstance(x, (int, float)) and x is not None:
                    out[f"{k}.{s}"] = float(x)
        elif isinstance(v, (int, float)) and v is not None:
            out[k] = float(v)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--tol", type=float, default=1e-6)
    ap.add_argument("--ffirst", action="store_true", help="also compare ffirst statistics")
    args = ap.parse_args()
    A = json.load(open(args.a, encoding="utf-8"))
    B = json.load(open(args.b, encoding="utf-8"))
    print(f"A = {A.get('impl')} accel={A.get('accel')}  ({args.a})")
    print(f"B = {B.get('impl')} accel={B.get('accel')}  ({args.b})")
    print(f"{'case':<12} {'coef':>9} {'se':>9} {'ffirst':>9}  counts  note")
    bad = 0
    for cid in A["cases"]:
        if cid not in B["cases"]:
            print(f"{cid:<12} missing in B")
            continue
        a, b = A["cases"][cid], B["cases"][cid]
        if not a.get("ok", True) or not b.get("ok", True):
            print(f"{cid:<12} error: A={a.get('error')} B={b.get('error')}")
            bad += 1
            continue
        rc, wc = maxrel(a.get("coef", {}), b.get("coef", {}))
        rs, ws = maxrel(a.get("se", {}), b.get("se", {}))
        rf, wf = (maxrel(ffirst_flat(a.get("ffirst")), ffirst_flat(b.get("ffirst")))
                  if args.ffirst else (0.0, ""))
        counts = all(a.get(k) == b.get(k) for k in ("N", "num_singletons", "num_separated"))
        fail = max(rc, rs, rf) > args.tol or not counts
        bad += fail
        note = ", ".join(x for x in (wc and f"coef:{wc}", ws and f"se:{ws}", wf and f"ff:{wf}") if x)
        if not counts:
            note += f" | N {a.get('N')} vs {b.get('N')}, sing {a.get('num_singletons')} vs " \
                    f"{b.get('num_singletons')}, sep {a.get('num_separated')} vs {b.get('num_separated')}"
        print(f"{cid:<12} {rc:9.2e} {rs:9.2e} {rf:9.2e}  {'ok' if counts else 'DIFF':<6}  "
              f"{'FAIL ' if fail else ''}{note}")
    print(f"{bad} case(s) outside tol {args.tol:g}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
