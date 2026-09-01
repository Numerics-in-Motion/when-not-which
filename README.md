# When, not which

A truck crosses a bridge. **The worst member is not one member.**

A Warren truss — **24 m span, 8 panels, 4 m deep**, a pin at one end and a
roller at the other. Every member is the same 100 mm square hollow section. A
two-axle lorry, **about 20 tonnes** (60 kN front, 140 kN rear, 4.0 m apart),
rolls across the deck.

At each of **241 positions the whole frame is solved again.** Nothing is
superposed and nothing is interpolated.

```
python -m pip install -r requirements.txt
python reproduce.py
```

Under a second.

## The answer

| the truck is | the worst member is |
|---|---|
| entering the span | an **end diagonal**, carrying the shear onto the bridge |
| anywhere in the middle | the **top chord** — and the worst chord **walks along the span, following the truck** |
| leaving | the **far diagonal** |

**It changes ten times on the way across**, over eight different members.

Peak **56 % of yield**, on a top chord, with the truck a little past midspan.

## And the reason the whole method exists

Park the truck at midspan, check the bridge there, and stop:

| | at midspan | at its own worst position | |
|---|---|---|---|
| one diagonal | **18 %** | **36 %** | **it doubles** |

**23 of the 31 members are under-read by more than five points of yield** that
way. That is what an influence line is for, and it is why a bridge is checked
with the load moving rather than parked.

## ⚠️ The check that does not trust any of this code

`m + 3 = 2j` — **31 + 3 = 34** — so the truss is statically determinate and the
support reactions are pure statics:

```
R_A = Σ Pᵢ (L − xᵢ) / L
```

**That expression does not read the element formulation, the member section,
or the stiffness assembly.** It would give the same answer for a truss made of
rubber. Change the section, change the material, change the element type — the
reactions must not move.

They do not: **2.4 × 10⁻¹⁴** at the worst of 241 positions.

This is worth saying plainly because the previous study on this rig had to
disclose the opposite. Its "independent" closed form read the same `E`, the
same `I` and the same fibre distance as the solver, so a wrong section property
would have moved both together — and once did, agreeing to 0.35 % while
carrying a 13.6 % error. **This check could catch that. That one could not.**

A second check needs no formula at all: under **one** axle the utilisation
curve must be symmetric about midspan, because the truss is. It is, to
**3.8 × 10⁻¹⁵**.

## What `reproduce.py` checks

| gate | |
|---|---|
| the shipped solver is byte-identical to the one that produced the reference | 7 files |
| the reactions match `Σ P(L−x)/L` at every position | worst 2.4e-14 |
| the truss is statically determinate, `m + 3 = 2j` | 31 + 3 = 34 |
| no position produced a mechanism | worst pivot ratio 1.5e-03 |
| under one axle the crossing is symmetric about midspan | worst 3.8e-15 |
| **the worst member changes, and both a diagonal and a chord take a turn** | 10 changes, 8 members |
| these frame elements are behaving as a truss | bending is 9.5 % of the peak |
| the bridge carries the truck everywhere on the span | peak 56.5 % |
| a midspan-only check materially under-reads at least three members | 23 of 31 |
| every number matches the frozen reference | 0.00e+00 |

**The sixth is pre-registered and could have killed the study.** If one member
had been worst from end to end there would be no *when* in the answer, and
this would have had to become a different video.

## On the number that is NOT quoted

The biggest **ratio** in the midspan-only comparison belongs to a different
member: it goes from 4.6 % to 21.4 %, which is a 78 % under-read. That number
is bigger and it is true, and it is not the headline, because **a member going
from 4.6 % to 21.4 % of yield is not a safety story** — quoting it would be
choosing the flattering measure. The claim here is the absolute one, 18
percentage points, and the count, 23 of 31. Neither can be improved by picking
a different member.

## What this is not

* **The elements carry bending as well as axial force**, which a textbook
  truss does not — the joints of a real bridge are not frictionless pins. The
  bending share of the peak stress is **measured rather than assumed away**:
  9.5 %.
* **Stress over yield is not a code check.** No load factors, no resistance
  factors, no buckling check on the compression members, no serviceability.
* **One truck, one lane, static.** No dynamic amplification, no braking force,
  no fatigue, no second vehicle, no wind.
* **No self-weight of the bridge itself** — this is the truck's effect alone,
  which is what an influence line isolates.
* One truss type, one span, one depth, one section.

## Files

```
reproduce.py                        the driver and the release gate
config.json                         span, panels, depth, section, truck, steel
bridge/                             the solver, byte-identical to the canonical one
  moving_load.py                    the crossing: axles, load placement, statics check
  geometry_bridges.py               the Warren truss generator and its supports
  geometry.py                       Structure, node de-duplication, load splitting
  fem.py                            frame assembly, constraints, member forces
  section_properties.py             the member section, with I and its fibre paired
  utils.py
reference/reference_canonical.json  frozen BEFORE the copy, with source hashes
results/results.json                written by reproduce.py
```
