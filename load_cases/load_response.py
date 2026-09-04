"""One bridge, two lorries, and the deflected shape of each.

WHAT THIS VIDEO HAS TO SHOW, AND WHY THIS MODULE EXISTS

The question being answered is the one asked most on this channel: what is
this computed with. A design review was blunt about how it can be answered
without falling into the band the channel's own numbers say is fatal:

    the promise that keeps it alive : "what happens to this bridge when the
                                      lorry gets heavier?"
    the promise that kills it       : "here is my simulation workflow"

So nothing here is a tour of software. It is the SAME bridge from the
published crossing study, at the SAME truck position, with one number
changed -- the axle loads -- solved twice. The software answers the question
by being the thing that made the second picture different from the first.

WHY A NEW MODULE RATHER THAN AN EDIT

`moving_load.solve_at` returns member stresses and throws the nodal
displacements away. The video needs the DEFLECTED SHAPE, which is the only
thing on screen a viewer can watch respond. That module ships inside a
published package and is frozen, so the displacements are taken here from the
same `fem.solve` call rather than by editing it.

WHAT MAY AND MAY NOT BE CLAIMED

The truss solve is LINEAR ELASTIC. There is no geometric nonlinearity, no
material failure, no buckling check on the compression members. So a heavier
lorry may honestly be shown to deflect the bridge further and stress it
harder, and may NOT be shown collapsing, yielding, or "failing" -- a design
review named that as the trap this subject invites.
"""
from __future__ import annotations

import copy
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import fem
import moving_load as ML

__all__ = ["with_axle_loads", "solve_shape", "midspan_deflection",
           "worst_at", "sweep_mass"]


def with_axle_loads(cfg, scale):
    """The same configuration with every axle load multiplied by `scale`.

    ONE NUMBER. The geometry, the section, the material, the axle spacing and
    the number of positions are untouched, so the two runs differ in exactly
    the quantity the video says they differ in.
    """
    c = copy.deepcopy(cfg)
    c["truck"]["axle_loads_n"] = [float(a) * float(scale)
                                  for a in cfg["truck"]["axle_loads_n"]]
    return c


def solve_shape(cfg, material, x_front):
    """One position, one load case: stresses AND the deflected node positions.

    Returns the solved `solve_at` dictionary plus `xy` -- the node positions
    in metres after displacement -- and the midspan deflection. The
    displacements come from the same `fem.solve` the published module calls,
    with the same assembled load vector, so the two are the same solve rather
    than two solves that ought to agree.
    """
    s, deck, _sec = ML.build_bridge(cfg)
    axles = ML.axle_positions(x_front, cfg)
    F = ML.load_vector(s, deck, axles)
    r = fem.solve(s, s.members, material["youngs_modulus_pa"], s.area,
                  s.inertia, F)
    u = np.asarray(r["u"], float).reshape(-1, 3)
    xy = np.asarray(s.nodes, float).copy()
    xy[:, 0] += u[:, 0]
    xy[:, 1] += u[:, 1]

    d = ML.solve_at(s, deck, cfg, material, x_front)
    d["xy"] = xy
    d["xy0"] = np.asarray(s.nodes, float).copy()
    d["u"] = u
    d["midspan_m"] = midspan_deflection(s, deck, u)
    return s, deck, d


def midspan_deflection(s, deck, u):
    """Downward movement of the deck node nearest midspan, in metres.

    POSITIVE MEANS DOWN, which is what a viewer expects a deflection to mean.
    The solver's y is up, so the sign is flipped here once and named, rather
    than being flipped at each place it is drawn -- a quantity whose sign
    depends on where you read it is how this channel has been wrong before.
    """
    nodes = np.asarray(s.nodes, float)
    span = float(nodes[:, 0].max() - nodes[:, 0].min())
    mid = 0.5 * (nodes[:, 0].max() + nodes[:, 0].min())
    cand = [i for i in deck] if len(deck) else list(range(len(nodes)))
    i = min(cand, key=lambda k: abs(nodes[k, 0] - mid))
    return float(-u[i, 1]), int(i), float(span)


def worst_at(d):
    """(member, utilisation) at this position, from the published solve."""
    return int(d["worst_member"]), float(d["worst_utilisation"])


def sweep_mass(cfg, material, scales, x_front):
    """The same position under several load scales.

    The video shows two. The sweep exists so the report can say whether the
    response is linear in the load -- which it must be, for a linear solve,
    and which is therefore a check on the assembly rather than a finding.
    """
    out = []
    for k in scales:
        c = with_axle_loads(cfg, k)
        _s, _deck, d = solve_shape(c, material, x_front)
        out.append(dict(scale=float(k),
                        total_n=float(sum(c["truck"]["axle_loads_n"])),
                        midspan_m=d["midspan_m"][0],
                        worst_member=int(d["worst_member"]),
                        worst_utilisation=float(d["worst_utilisation"])))
    return out
