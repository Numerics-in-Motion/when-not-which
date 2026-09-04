"""The same bridge, the same lorry position, one number changed.

    python load_cases/reproduce_load_cases.py
    python load_cases/reproduce_load_cases.py --no-check

WHAT THIS ADDS TO THE CROSSING STUDY BESIDE IT

The crossing study rolls one lorry across and asks which member is worst and
when. This asks a smaller question with the same solver: change the axle
loads in `config.json`, leave everything else alone, and see what the bridge
does.

    200 kN (20.4 t)   midspan 25.16 mm   worst member 52.4 % of yield
    300 kN (30.6 t)   midspan 37.74 mm   worst member 78.6 % of yield

WHY NOT DOUBLE IT

Doubling the lorry reaches 104.7 % of yield. This solve is LINEAR ELASTIC --
no geometric nonlinearity, no material failure, no buckling check on the
compression members -- so past yield it is describing nothing. The
multiplier is 1.5 because that is what the model supports.

**Nothing here yields and nothing collapses.**

THE ONE CHECK THAT DOES NOT READ THE SOLVER

The truss is statically determinate, so the reactions are

    R_A = sum P_i (L - x_i) / L

which does not know the element formulation, the section, or how the
stiffness was assembled. It would give the same answer for a truss made of
rubber. It matches to about 1e-14. Everything else in here -- including the
exact proportionality between load and deflection -- is the solver agreeing
with itself: a check that the load vector scaled and nothing else did, not
evidence about bridges.

CHANGE SOMETHING

`config.json`, `truck.axle_loads_n`. Two numbers. Then run this again --
and note that past about 1.9x the answer stops meaning anything, for the
reason above.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(PKG, "bridge"))
sys.path.insert(0, HERE)

import load_response as LR                                     # noqa: E402

REF = os.path.join(HERE, "reference_load_cases.json")
CONFIG = os.path.join(PKG, "config.json")


def sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def gate(label, ok, detail=""):
    print("   [%s] %s%s" % ("PASS" if ok else "FAIL", label,
                            ("   %s" % detail) if detail else ""))
    return bool(ok)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-check", action="store_true")
    a = ap.parse_args(argv)

    with open(REF, encoding="utf-8") as fh:
        ref = json.load(fh)
    with open(CONFIG, encoding="utf-8") as fh:
        cfg = json.load(fh)
    mat = cfg["material"]
    x = float(ref["meta"]["x_front_m"])

    print(__doc__.splitlines()[0])
    print()
    t0 = time.time()

    out = {}
    for k in ref["meta"]["scales"]:
        c = LR.with_axle_loads(cfg, k)
        t1 = time.perf_counter()
        s, deck, d = LR.solve_shape(c, mat, x)
        dt = time.perf_counter() - t1
        mid, mid_node, span = d["midspan_m"]
        out["%.2f" % k] = dict(
            total_n=float(sum(c["truck"]["axle_loads_n"])),
            midspan_m=float(mid), worst_member=int(d["worst_member"]),
            worst_utilisation=float(d["worst_utilisation"]),
            reaction_a=float(d["reaction_a"]),
            reaction_a_exact=float(d["reaction_a_exact"]),
            solve_seconds=float(dt))
        print("   %6.0f kN   midspan %6.2f mm   member %2d at %5.1f %% of "
              "yield   (%.0f ms)"
              % (out["%.2f" % k]["total_n"] / 1e3, 1000 * mid,
                 d["worst_member"], 100 * d["worst_utilisation"], 1000 * dt))

    lo, hi = out["1.00"], out["1.50"]
    res = dict(cases=out, x_front_m=x,
               ratio_deflection=hi["midspan_m"] / lo["midspan_m"],
               meta=ref["meta"])
    d = os.path.join(HERE, "results")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "load_cases.json"), "w", encoding="utf-8",
              newline=chr(10)) as fh:
        json.dump(res, fh, indent=1)

    if a.no_check:
        print("\n   (release gate skipped)")
        return 0

    print()
    ok = True
    for mod, want in sorted(ref["source_sha256"].items()):
        got = sha(os.path.join(HERE, mod + ".py"))
        if got != want:
            print("       %s.py differs" % mod)
            ok = False
    ok = gate("the shipped module is byte-identical to the one that produced "
              "the reference", ok, "%d file" % len(ref["source_sha256"]))
    ok &= gate("config.json is the one the reference was frozen with",
               sha(CONFIG) == ref["config_sha256"])

    for k in ("1.00", "1.50"):
        w = ref["cases"][k]
        ok &= gate("%s: midspan reproduces" % k,
                   abs(out[k]["midspan_m"] - w["midspan_m"])
                   / w["midspan_m"] < 1e-9,
                   "%.6f mm" % (1000 * out[k]["midspan_m"]))

    ok &= gate("both cases stay below yield -- the reason the multiplier is "
               "1.5 and not 2",
               all(v["worst_utilisation"] < 0.95 for v in out.values()),
               " ".join("%.1f %%" % (100 * v["worst_utilisation"])
                        for v in out.values()))

    err = max(abs(v["reaction_a"] - v["reaction_a_exact"])
              / abs(v["reaction_a_exact"]) for v in out.values())
    ok &= gate("the reactions match pure statics -- the only check here that "
               "does not read the solver", err < 1e-9, "%.2e" % err)

    ratio = hi["midspan_m"] / lo["midspan_m"]
    ok &= gate("the response is proportional (a property of a linear solve, "
               "not a finding)", abs(ratio - 1.5) < 1e-9, "%.9f" % ratio)

    # and the one that would catch someone quietly turning the load up
    c2 = copy.deepcopy(cfg)
    c2["truck"]["axle_loads_n"] = [v * 2.0
                                   for v in cfg["truck"]["axle_loads_n"]]
    _s, _deck, d2 = LR.solve_shape(c2, mat, x)
    ok &= gate("and doubling it really would leave the model behind",
               d2["worst_utilisation"] > 1.0,
               "%.1f %% of yield at x2" % (100 * d2["worst_utilisation"]))

    print("\n   %s in %.1f s" % ("ALL PASS" if ok else "FAILED",
                                 time.time() - t0))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
