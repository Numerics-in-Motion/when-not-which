"""Geometry generation for the three competing bridge structures.

All three bridges share the SAME span, the SAME height envelope (used only
for a common view-box / bounding height so panels compare visually), the
SAME support span, the SAME midspan load point and the SAME uniform-area
material-volume matching rule as the shapes video (spec 5.2/5.4). They
differ only in how they carry the load from midspan to the two supports:

    beam  : a single straight line, simply supported -> carries load almost
            entirely by bending (the weak case, by design).
    truss : a Warren truss (top/bottom chord + zigzag diagonals, no
            verticals) -> carries load mostly as chord/diagonal axial force,
            continuing the "triangulation wins" story from video 1.
    arch  : a two-hinged parabolic arch rib -> carries load mostly as
            compressive axial thrust, resisted by BOTH supports being full
            pins (see `roller_is_pin` on `Structure` / `fem.constrained_dofs`).
            This is the physically correct support for an arch: a roller at
            one end would not develop the horizontal thrust that makes an
            arch strong, so the two structures would not be comparing "arch
            behaviour" at all. Every structure still shares the same pin
            position, span and load line -- only the *type* of the second
            support differs, and only for the arch.

The reused pieces from geometry.py: `Structure`, `_nearest_node`,
`_split_load_at_x` (load resultant placement) and `assign_matched_area`
(uniform member area from total length + target volume) are identical to
the shapes video -- nothing about the FEM engine, failure model or the
material-volume rule is reimplemented here.
"""
from __future__ import annotations

import numpy as np

from geometry import (Structure, _nearest_node, _split_load_at_x,
                      assign_matched_area)


def generate_beam(span: float, view_height: float, n_panels: int = 10) -> tuple:
    """Straight simply-supported beam, subdivided into `n_panels` segments.
    The whole line is the "deck" -- every node is a load candidate."""
    xs = np.linspace(0.0, span, n_panels + 1)
    nodes = np.array([(x, 0.0) for x in xs])
    members = [(i, i + 1) for i in range(n_panels)]
    s = Structure("beam", nodes, members, span, view_height)
    return s, np.arange(s.n_nodes)


def generate_truss_warren(span: float, truss_height: float, n_panels: int,
                          view_height: float) -> tuple:
    """Warren truss: bottom + top chord with zigzag diagonals, no verticals.

    `n_panels` (bottom-chord panel count) is forced to even so an exact
    midspan node exists on the bottom chord, which is treated as the deck.
    """
    if n_panels % 2 != 0:
        n_panels += 1
    bx = np.linspace(0.0, span, n_panels + 1)
    bottom = [(x, 0.0) for x in bx]
    panel = span / n_panels
    tx = (np.arange(n_panels) + 0.5) * panel
    top = [(x, truss_height) for x in tx]
    nodes = np.array(bottom + top)
    nb = n_panels + 1

    members = []
    for i in range(n_panels):
        members.append((i, i + 1))                  # bottom chord
    for i in range(n_panels - 1):
        members.append((nb + i, nb + i + 1))         # top chord
    for i in range(n_panels):
        members.append((i, nb + i))                  # diagonal, rising
        members.append((nb + i, i + 1))               # diagonal, falling

    s = Structure("truss", nodes, members, span, view_height)
    return s, np.arange(n_panels + 1)  # bottom chord = deck


def generate_arch(span: float, rise: float, n_segments: int,
                  view_height: float) -> tuple:
    """Two-hinged parabolic arch rib, y = rise * (1 - ((2x/L)-1)^2),
    approximated by `n_segments` straight frame elements. The whole rib is
    the deck (load applied directly to the arch, as in video 1's "load on
    the structure being analyzed" convention)."""
    xs = np.linspace(0.0, span, n_segments + 1)
    ys = rise * (1.0 - ((2.0 * xs / span) - 1.0) ** 2)
    nodes = np.array(list(zip(xs, ys)))
    members = [(i, i + 1) for i in range(n_segments)]
    s = Structure("arch", nodes, members, span, view_height)
    s.roller_is_pin = True  # two-hinged: both supports resist horizontal thrust
    return s, np.arange(s.n_nodes)


def set_boundary_conditions_bridge(s: Structure, load_candidate_idx) -> None:
    """Pin at (0,0), roller (or second pin, for the arch) at (span,0); load
    resultant placed at midspan among the given deck-candidate nodes using
    the same `_split_load_at_x` helper as the shapes video."""
    L = s.width
    s.pin_node = _nearest_node(s.nodes, (0.0, 0.0))
    s.roller_node = _nearest_node(s.nodes, (L, 0.0))
    s.load_nodes, s.load_weights = _split_load_at_x(
        s.nodes, np.asarray(load_candidate_idx), L / 2.0)


def build_all_bridges(config: dict) -> list:
    """Generate beam/truss/arch, assign BCs and matched member areas."""
    g = config["geometry"]
    span = float(g["span_m"])
    h_max = float(g["max_height_m"])
    target_v = float(g["target_material_volume_m3"])

    n_beam = int(g.get("beam_n_panels", 10))
    n_truss = int(g.get("truss_n_panels", 8))
    n_arch = int(g.get("arch_n_segments", 24))
    truss_height = float(g.get("truss_height_ratio", 1.0)) * h_max
    arch_rise = float(g.get("arch_rise_ratio", 0.6)) * h_max

    generated = [
        generate_beam(span, h_max, n_beam),
        generate_truss_warren(span, truss_height, n_truss, h_max),
        generate_arch(span, arch_rise, n_arch, h_max),
    ]
    structures = []
    for s, load_candidates in generated:
        set_boundary_conditions_bridge(s, load_candidates)
        assign_matched_area(s, target_v)
        structures.append(s)
    return structures
