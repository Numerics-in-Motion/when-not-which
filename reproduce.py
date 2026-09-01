"""A truck crosses a bridge. Which member is worst, and when?

    python -m pip install -r requirements.txt
    python reproduce.py                 # about a minute
    python reproduce.py --no-check      # run without the release gate

A Warren truss, 24 m span, 8 panels, 4 m deep, simply supported: a pin at one
end and a roller at the other. Every member is the same 100 mm square hollow
section. A two-axle lorry -- 60 kN on the front axle, 140 kN on the rear,
4.0 m apart, about 20 tonnes -- rolls across the deck.

At each of 241 positions the whole frame is solved again. Nothing is
superposed and nothing is interpolated.

THE ANSWER

The worst member is not one member. It changes ten times on the way across:

    entering        an end diagonal, carrying the shear onto the span
    the middle      the top chord -- and the worst chord WALKS along the
                    span, following the truck
    leaving         the far diagonal

Peak 56 % of yield, on a top chord, with the truck a little past midspan.

AND THE POINT OF THE WHOLE METHOD

Park the truck at midspan, check the bridge there, and stop -- and one
diagonal reads 18 % when it actually reaches 36 %. It DOUBLES. Twenty-three
of the thirty-one members are under-read by more than five points of yield
that way.

That is what an influence line is for, and it is why a bridge is checked with
the load moving rather than parked.

WHY THE ANSWER CAN BE CHECKED WITHOUT TRUSTING ANY OF THIS CODE

`m + 3 = 2j`, so the truss is statically determinate and the support
reactions are pure statics:

    R_A = sum P_i (L - x_i) / L

That expression does not read the element formulation, the member section, or
the stiffness assembly. It would be the same for a truss made of rubber.
Change the section, change the material, change the element type -- the
reactions must not move. They do not: the solver matches to about 1e-14 at
every one of the 241 positions.

A second check needs no formula at all. Under ONE axle the utilisation curve
must be symmetric about midspan, because the truss is.

WHAT THIS IS NOT

The elements carry bending as well as axial force, which a textbook truss
does not -- the joints of a real bridge are not frictionless pins. The
bending share of the peak stress is measured rather than assumed away; it
comes out at 9.5 %.

One truck, one lane, static. No dynamic amplification, no braking force, no
fatigue, no second vehicle, no wind, and no self-weight of the bridge itself.
Stress over yield is not a code check: there are no load factors, no
resistance factors and no buckling check on the compression members here.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "bridge"))

import moving_load as ml                                       # noqa: E402

REF = os.path.join(HERE, "reference", "reference_canonical.json")


def sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def gate(label, ok, detail=""):
    print("   [%s] %s%s" % ("PASS" if ok else "FAIL", label,
                            ("   %s" % detail) if detail else ""))
    return bool(ok)


def check_sources(ref):
    ok = True
    for mod, want in sorted(ref["source_sha256"].items()):
        got = sha(os.path.join(HERE, "bridge", mod + ".py"))
        if got != want:
            print("       %s.py differs: %s != %s" % (mod, got[:12], want[:12]))
            ok = False
    ok &= sha(os.path.join(HERE, "config.json")) == ref["config_sha256"]
    return gate("the shipped solver is byte-identical to the one that "
                "produced the reference", ok,
                "%d files" % (len(ref["source_sha256"]) + 1))


def runs(hist):
    out = []
    for h in hist:
        if out and out[-1][0] == h["worst_member"]:
            out[-1][2] = h["x_front"]
        else:
            out.append([h["worst_member"], h["x_front"], h["x_front"]])
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-check", action="store_true",
                    help="run without comparing to the frozen reference")
    a = ap.parse_args(argv)
    t0 = time.time()

    with open(REF, encoding="utf-8") as f:
        ref = json.load(f)
    with open(os.path.join(HERE, "config.json"), encoding="utf-8") as f:
        cfg = json.load(f)
    mat = cfg["material"]
    span = float(cfg["geometry"]["span_m"])

    print("A truck crosses a bridge. Which member is worst, and when?")
    print("   %.0f m span, %d panels, %.0f m deep, %s"
          % (span, cfg["geometry"]["n_panels"],
             cfg["geometry"]["truss_height_m"], cfg["truck"]["name"]))
    print("   axles %s m apart, %s N\n"
          % (cfg["truck"]["axle_offsets_m"], cfg["truck"]["axle_loads_n"]))

    s, deck, sec, hist = ml.crossing(cfg, mat)
    kinds = [ml.member_kind(s, k, deck) for k in range(len(s.members))]

    seq = runs(hist)
    print("   the worst member, position by position")
    for member, x0, x1 in seq:
        print("      %5.1f - %5.1f m   member %-3d  %s"
              % (x0, x1, member, kinds[member]))

    peak = max(hist, key=lambda h: h["worst_utilisation"])
    print("\n   peak %.1f %% of yield at x = %.2f m, member %d (%s)"
          % (100 * peak["worst_utilisation"], peak["x_front"],
             peak["worst_member"], kinds[peak["worst_member"]]))

    # the envelope, and what a midspan-only check misses
    env = np.zeros(len(s.members))
    mid_i = min(range(len(hist)),
                key=lambda i: abs(hist[i]["x_front"] - span / 2.0))
    for h in hist:
        for r in h["members"]:
            env[r["member"]] = max(env[r["member"]], r["utilisation"])
    mid = np.array([r["utilisation"] for r in
                    sorted(hist[mid_i]["members"], key=lambda r: r["member"])])
    gapv = env - mid
    k = int(np.argmax(gapv))
    print("\n   midspan-only check")
    print("      member %d (%s): %.1f %% at midspan, %.1f %% at its own worst"
          % (k, kinds[k], 100 * mid[k], 100 * env[k]))
    print("      that is %.1f points, a factor of %.2f"
          % (100 * gapv[k], env[k] / max(mid[k], 1e-12)))
    print("      %d of the %d members are under-read by more than 5 points"
          % (int(np.sum(gapv > 0.05)), len(env)))

    print()
    if a.no_check:
        print("   (release gate skipped)")
        return 0

    ok = check_sources(ref)

    # -- statics. The check that does not read the solver's own inputs. -----
    rel = 0.0
    for h in hist:
        for got, want in ((h["reaction_a"], h["reaction_a_exact"]),
                          (h["reaction_b"], h["reaction_b_exact"])):
            if want > 1.0:
                rel = max(rel, abs(got - want) / want)
    ok &= gate("the reactions match sum P(L-x)/L at every position", rel < 1e-9,
               "worst %.2e" % rel)

    m, j = len(s.members), int(s.n_nodes)
    ok &= gate("the truss is statically determinate, m + 3 = 2j",
               m + 3 == 2 * j, "%d + 3 = %d" % (m, 2 * j))
    ok &= gate("no position produced a mechanism",
               not any(h["singular"] for h in hist),
               "worst pivot ratio %.2e"
               % min(h["min_pivot_ratio"] for h in hist))

    # -- symmetry. Needs no formula at all. ---------------------------------
    one = dict(cfg)
    one["truck"] = {"name": "single axle", "axle_offsets_m": [0.0],
                    "axle_loads_n": [float(sum(cfg["truck"]["axle_loads_n"]))]}
    _s, _d, _sec, h1 = ml.crossing(one, mat, n_steps=121)
    on = [h for h in h1 if h["x_front"] <= span + 1e-9]
    u = np.array([h["worst_utilisation"] for h in on])
    asym = float(np.max(np.abs(u - u[::-1])) / np.max(u))
    ok &= gate("under one axle the crossing is symmetric about midspan",
               asym < 1e-6, "worst %.2e" % asym)

    # -- the pre-registered one ---------------------------------------------
    ks = [kinds[r[0]] for r in seq]
    ok &= gate("the worst member changes, and both a diagonal and a chord "
               "take a turn",
               len(seq) - 1 >= 2 and "diagonal" in ks
               and any("chord" in x for x in ks),
               "%d changes, %d distinct members"
               % (len(seq) - 1, len({r[0] for r in seq})))

    share = max(h["bending_share"] for h in hist)
    ok &= gate("these frame elements are behaving as a truss", share < 0.20,
               "bending is %.1f %% of the peak stress" % (100 * share))

    ok &= gate("the bridge carries the truck everywhere on the span",
               peak["worst_utilisation"] < 1.0,
               "peak %.1f %%" % (100 * peak["worst_utilisation"]))

    ok &= gate("a midspan-only check materially under-reads at least three "
               "members", gapv[k] > 0.05 and int(np.sum(gapv > 0.05)) >= 3,
               "%d of %d, worst %.1f points"
               % (int(np.sum(gapv > 0.05)), len(env), 100 * gapv[k]))

    # -- and against the frozen numbers -------------------------------------
    dev = max(abs(env[i] - ref["envelope"][i]) for i in range(len(env)))
    ok &= gate("every number matches the frozen reference", dev < 1e-9,
               "worst %.2e" % dev)

    out = os.path.join(HERE, "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "results.json"), "w", encoding="utf-8",
              newline="\n") as f:
        json.dump({"span_m": span,
                   "members": [{"member": i, "kind": kinds[i],
                                "envelope": float(env[i]),
                                "midspan": float(mid[i])}
                               for i in range(len(env))],
                   "migration": [{"member": r[0], "kind": kinds[r[0]],
                                  "from_m": r[1], "to_m": r[2]} for r in seq],
                   "peak": {"utilisation": peak["worst_utilisation"],
                            "x_front_m": peak["x_front"],
                            "member": peak["worst_member"]},
                   "midspan_only": {"member": k, "kind": kinds[k],
                                    "midspan": float(mid[k]),
                                    "envelope": float(env[k]),
                                    "points": float(gapv[k]),
                                    "n_underread_5_points":
                                        int(np.sum(gapv > 0.05))}},
                  f, indent=1)

    print("\n   %s in %.1f s" % ("ALL GATES PASS" if ok else "GATES FAILED",
                                 time.time() - t0))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
