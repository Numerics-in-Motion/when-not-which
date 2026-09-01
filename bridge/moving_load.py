"""A truck crosses a bridge. Which member becomes critical, and when?

WHY THIS EXISTS

Every study on this rig for the last four videos asked the same shape of
question -- "which of these four is strongest" -- of the same object, a
column. This asks a different KIND of question of a different object: not
*which* member is worst, but *when* it is, and the answer is that the worst
member is not the same member for the whole crossing.

    near the abutments   the end diagonals carry the shear
    near midspan         the chords carry the moment

So the critical member migrates as the truck rolls, comes back, and the
structure has no single weakest link -- it has a weakest link per position.
That is an influence line, and it is why bridge engineers move the load
instead of placing it in the middle and stopping.

WHAT IS SOLVED

A Warren truss, simply supported, span and panels from config. A two-axle
truck rolls across the deck. At every position the whole frame is re-solved
and every member's combined stress is recorded. Nothing is superposed or
interpolated: each position is its own solve.

WHY THE ANSWER CAN BE CHECKED EXACTLY

The truss is **statically determinate** (m + 3 = 2j), so:

  * the support reactions are pure statics -- `R_A = sum P_i (L - x_i) / L` --
    and do not depend on the element formulation at all. The solver matches
    that to about 1e-14, which is a real external check rather than the rig
    agreeing with itself;
  * a symmetric structure under a single load must give an influence line
    symmetric about midspan, which is a second check that needs no formula.

WHAT THIS IS NOT

The elements are FRAME elements, so they carry bending as well as axial
force; a textbook truss carries only axial force. That is deliberate -- the
joints of a real bridge are not frictionless pins -- and the bending share of
the peak stress is measured and reported rather than assumed small. It comes
out at a few percent.

One truck, one lane, static. No dynamic amplification, no braking force, no
fatigue, no second vehicle, no wind, and no self-weight of the bridge.
"""
from __future__ import annotations

import numpy as np

import fem
import section_properties as sp
from geometry import _split_load_at_x
from geometry_bridges import generate_truss_warren, set_boundary_conditions_bridge


def build_bridge(cfg):
    """The truss, its supports, and one member section shared by all members.

    The section is a square hollow tube solved by `section_properties`, not a
    hand-written A and I: the pairing of the second moment with the fibre
    distance it belongs to is exactly the thing this channel got wrong once,
    and there is no reason to write it out a second time by hand.
    """
    g = cfg["geometry"]
    span = float(g["span_m"])
    s, deck = generate_truss_warren(span, float(g["truss_height_m"]),
                                    int(g["n_panels"]), float(g["truss_height_m"]))
    set_boundary_conditions_bridge(s, deck)
    sec = sp.box_properties(float(g["member_area_m2"]),
                            float(g["member_outer_m"]))
    s.area = float(sec["area"])
    s.inertia = float(sec["I_min"])
    s.fibre_c = sp.fibre_for(sec, sec["I_min"])
    return s, np.asarray(deck), sec


def axle_positions(x_front, cfg):
    """(position, load) for each axle, keeping only the ones on the deck."""
    t = cfg["truck"]
    out = []
    for offset, load in zip(t["axle_offsets_m"], t["axle_loads_n"]):
        x = x_front - float(offset)
        if -1e-9 <= x <= float(cfg["geometry"]["span_m"]) + 1e-9:
            out.append((float(np.clip(x, 0.0, cfg["geometry"]["span_m"])),
                        float(load)))
    return out


def load_vector(s, deck, axles):
    """Axle loads shared between the two deck nodes either side, by lever rule.

    `_split_load_at_x` is the same helper the published bridge study used to
    place its midspan load, so a truck between panel points is placed the way
    this rig has always placed a load.
    """
    F = np.zeros(3 * s.n_nodes)
    for x, P in axles:
        nodes, w = _split_load_at_x(s.nodes, deck, x)
        for n, wt in zip(nodes, w):
            F[3 * n + 1] -= P * wt
    return F


def exact_reactions(axles, span):
    """Statics. Independent of the element formulation, and therefore a check.

    A determinate, simply supported span carries `P (L - x) / L` to the near
    support and `P x / L` to the far one, whatever the members are made of or
    how they are connected.
    """
    ra = sum(P * (span - x) / span for x, P in axles)
    rb = sum(P * x / span for x, P in axles)
    return ra, rb


def solve_at(s, deck, cfg, material, x_front):
    """One truck position: every member's stress, and the worst one."""
    axles = axle_positions(x_front, cfg)
    F = load_vector(s, deck, axles)
    r = fem.solve(s, s.members, material["youngs_modulus_pa"], s.area,
                  s.inertia, F)
    fy = material["failure_stress_pa"]
    rows = []
    for (m, i, j, L, axial, mi, mj) in r["member_results"]:
        s_ax = abs(axial) / s.area
        s_bd = max(abs(mi), abs(mj)) * s.fibre_c / s.inertia
        rows.append(dict(member=m, node_i=i, node_j=j, length=L,
                         axial_n=axial, axial_pa=s_ax, bending_pa=s_bd,
                         combined_pa=s_ax + s_bd,
                         utilisation=(s_ax + s_bd) / fy))
    worst = max(rows, key=lambda d: d["combined_pa"])
    ra_exact, rb_exact = exact_reactions(axles, cfg["geometry"]["span_m"])
    return dict(x_front=x_front, axles=axles, members=rows,
                worst_member=worst["member"],
                worst_utilisation=worst["utilisation"],
                bending_share=(worst["bending_pa"] / worst["combined_pa"]
                               if worst["combined_pa"] > 0 else 0.0),
                reaction_a=float(r["reactions"][3 * s.pin_node + 1]),
                reaction_b=float(r["reactions"][3 * s.roller_node + 1]),
                reaction_a_exact=ra_exact, reaction_b_exact=rb_exact,
                singular=bool(r["singular"]),
                min_pivot_ratio=float(r["min_pivot_ratio"]))


def crossing(cfg, material, n_steps=None):
    """The whole crossing, front axle from the near abutment to clear of it."""
    s, deck, sec = build_bridge(cfg)
    span = float(cfg["geometry"]["span_m"])
    tail = max(cfg["truck"]["axle_offsets_m"])
    n = int(n_steps or cfg["solver"]["n_positions"])
    xs = np.linspace(0.0, span + tail, n)
    return s, deck, sec, [solve_at(s, deck, cfg, material, float(x)) for x in xs]


def member_kind(s, m, deck):
    """bottom chord (the deck) / top chord / diagonal, from the geometry.

    Derived from the node coordinates rather than from the order the
    generator happened to append them in, so a change of generator cannot
    silently relabel the answer.
    """
    i, j = s.members[m]
    yi, yj = s.nodes[i, 1], s.nodes[j, 1]
    if abs(yi) < 1e-9 and abs(yj) < 1e-9:
        return "deck"
    if abs(yi - yj) < 1e-9:
        return "top chord"
    return "diagonal"
