# L-Shape Re-entrant Corner Meshing Issue

Tracking the investigation and fix for non-straight boundary edges at the
re-entrant (concave) corner of an L-shaped domain meshed with `pyPolyMesher`.

- **Status:** ✅ FIXED — paper's method implemented as `reentrant_corner_seeds()`
- **Domain:** L-shape = 2×2 square minus top-right 1×1 square; re-entrant corner at `(1, 1)`
- **Regression tests:** `tests/test_reentrant_corner.py`, `tests/test_mesh_quality.py`

---

## 1. Problem Statement

When meshing the L-shaped domain, the boundary edges around the inner
(re-entrant) corner `(1, 1)` are **not straight**, even though a mesh node is
pinned exactly at the corner via `PFix`. The boundary makes a small notch/kink
right next to the apex.

---

## 2. Investigation Log

### Step 1 — Confirm the corner node exists
- A node sits **exactly** at `(1, 1)` (distance `0.0`) once `[1,1]` is added to `PFix`.
- So the apex itself is fine; the defect is elsewhere.

### Step 2 — Inspect boundary edges out of the apex (`dPolygon` SDF)
Tracing edges from corner node `(1,1)`:
- **Vertical wall:** `(1,1) → (1.000,1.213) → (1.000,1.434)` ✅ straight
- **Horizontal wall:** `(1,1) → (1.118, 0.993) → (1.126, 1.000)` ❌
  - Node `(1.118, 0.993)` dips **0.007 inside** the domain → notch/kink.

### Step 3 — Rule out the SDF
- Switched SDF from `dPolygon` to `dDiff(dRectangle, dRectangle)`.
- The wobble **moved to the other wall** instead of disappearing; one node was
  even pushed slightly **outside** (into the notch).
- **Conclusion:** the distance-function choice is *not* the cause.

### Step 4 — Characterize as stochastic, corner-intrinsic
- Metric: max deviation of near-corner boundary nodes from the notch walls.
- Across 12 random seeds: **7/12 produced a kink** (deviation > 0.005), worst `0.13`.
- **Root cause:** A PolyMesher boundary edge is straight only when it is the
  Voronoi edge between a seed and its mirror image across that wall. Reflection
  guarantees this *away* from corners, but the boundary node **adjacent** to a
  concave apex is pinned by free random seeds whose symmetry across the wall is
  not guaranteed. Pinning only the apex does not constrain its neighbors.

---

## 2b. Reference Paper (authoritative diagnosis & fix)

Source: Talischi, Paulino, Pereira, Menezes, *"PolyMesher: a general-purpose mesh
generator for polygonal elements written in Matlab"*, Struct Multidisc Optim (2012)
— `SMO_12_PolyMesher.pdf`, **Fig. 15** and surrounding text. This is exactly our issue.

**The paper confirms the root cause:**
> "Since the seeds near the non-convex corner are placed independently (during
> the Lloyd's iterations), their reflections may interfere with each other."

**Fig. 15 caption:**
> "In the left figure, the reflection of point y1 has interfered with the
> approximation of the horizontal boundary by seeds z and y2. In the right
> figure, seeds y1, y2 and y3 are fixed in such a way that they all have y as
> the reflection."

**Prescribed fix (canonical):**
> "One possible solution is to place some seeds at an equal distance from the
> corner and fix them during the Lloyd's iteration. This process can be made
> systematic by choosing, as the fixed seeds, the reflections of an
> appropriately placed seed outside such a corner (e.g. point y in Fig. 15)."

**Related (Fig. 7, convex corners):**
> "To accurately capture a corner, nearby seeds need to be reflected about both
> boundary segments incident on that corner."

**Interpretation for our L-shape (corner at `(1,1)`, notch above-right):**
1. Pick one external point `y` placed *outside* the re-entrant corner (inside the notch),
   e.g. `y = (1+δ, 1+δ)`.
2. Reflect `y` about each of the two notch walls (`x=1` and `y=1`) to get the
   fixed seeds `y1 = (1-δ, 1+δ)` and `y2 = (1+δ, 1-δ)`; optionally a third
   `y3` so all share the same reflection `y`.
3. Add `y1, y2, y3` to `PFix` so they stay put during Lloyd's iterations.
4. Because the fixed seeds are equidistant from the corner and share the common
   external reflection `y`, their reflections no longer interfere and both notch
   walls are approximated straight up to the apex.

Note: this is listed in the paper as an **extension deliberately left out of the
released Matlab code** "for the sake of simplicity and clarity" — so the stock
PolyMesher (and this Python port) does **not** implement it. That is why pinning
only the apex is insufficient.

---

## 3. Fixes Tested

| Approach | Result |
|---|---|
| Pin apex only | 7/12 kinked, worst `0.13` |
| Pin seeds *on* the wall lines (`x=1`, `y=1`) | ❌ degenerate generators, worse (`0.19`) |
| Pin seeds *just inside* each notch wall | worst `0.13 → 0.03`, but still frequent small kinks (11/12) |
| **Post-process: snap near-wall boundary nodes onto the wall** | ✅ **7/12 → 1/12 kinked** |

### Snapping detail
- Project each boundary node within ~½ element-size of a wall exactly onto that wall.
- `tol = 0.25 * sqrt(area/NElem)` → kinks 7/12 → 1/12.
- Residual failure: a rare seed where an **entire element** lands in the notch —
  a topological defect, not a small drift; snapping cannot rescue it.

---

## 4. Implementation (DONE)

- [x] **Primary fix (paper's canonical method, Fig. 15):** added
  `reentrant_corner_seeds(corner, e1, e2, delta)` in
  `src/pyPolyMesher/pyPolyMesher.py`, exported from `pyPolyMesher`.
  - Places external point `y = corner + delta*(e1+e2)` in the notch.
  - Returns the three fixed seeds = reflections of `y` about each incident
    wall and through the corner point (`y1, y2, y3`), all equidistant from the
    corner and all reflecting back onto `y`.
  - Append the result to `Domain.PFix`; seeds persist through Lloyd's iterations.
- [x] **Recommended usage:** keep the apex `[1,1]` pinned **and** add the three
  re-entrant seeds with `delta ~ 1.6 * element_size`
  (`element_size = sqrt(area/NElem)`). Encoded in `tests/test_reentrant_corner.py`.
  - NOTE: omitting the apex seed makes it worse — pin the apex too.
  - `delta` too small (<1*h) lets free seeds crowd in; the 1.4-2.0*h band
    captures the apex exactly.
- [x] Verified across 12 random seeds (see section 5).
- [ ] **Optional polish:** boundary-snapping post-process to remove the residual
  sub-element wobble (~0.01-0.02). Not required now that the catastrophic case
  is gone.
- [ ] **Secondary:** smaller elements — absolute kink scales with element size.
- [ ] Add a permanent regression test asserting near-corner boundary nodes lie
  on the walls and the apex is captured.

### API
```python
from pyPolyMesher import reentrant_corner_seeds
import numpy as np

h = np.sqrt(area / NElem)              # ~ element size
PFix = [[1, 1], [0, 0], [2, 0], [2, 1], [1, 2], [0, 2]]   # apex + convex corners
PFix += reentrant_corner_seeds(corner=[1, 1], e1=[1, 0], e2=[0, 1], delta=1.6 * h)
```

---

## 5. Key Metrics Reference

Max deviation of near-corner boundary nodes from the notch walls, over random seeds.

| Config | mean dev | worst dev | apex captured |
|---|---|---|---|
| apex only (before) | 0.01622 | 0.12966 | always (pinned) |
| inside-wall seeds (ad-hoc) | 0.01394 | 0.03339 | — |
| **paper fix, apex + y1y2y3, d=1.6*h (after)** | **0.0097** | **0.0197** | **12/12 exact** |
| paper fix, d=2.0*h | 0.0056 | 0.0187 | 12/12 exact |
| apex + snap (tol 0.25*h) | 0.01081 | 0.12966 | 1/12 cleared |

Catastrophic single-seed case (seed 4): `0.1297 -> 0.0112` (12x better), apex gap `0.0366 -> 0.0000`.

`h = sqrt(area/NElem) = sqrt(3/100) ~ 0.173`

---

## 6. Files

- `src/pyPolyMesher/pyPolyMesher.py` — `reentrant_corner_seeds()` helper (the fix)
- `src/pyPolyMesher/__init__.py` — exports `reentrant_corner_seeds`
- `tests/test_reentrant_corner.py` — regression tests for the re-entrant seeds
- `tests/test_mesh_quality.py` — randomized L-shape robustness tests (the former
  `mesh_lshape.py` / `validate_lshape.py` repro scripts, retired into tests)
- `docs/LSHAPE_CORNER_ISSUE.md` — this tracking document

Note: package is now installed editable (`pip install -e .`) so `src/` is live.

---

## 7. Critical Bug: nodes flung far outside the domain (FIXED)

**Symptom:** final meshes contained element vertices far outside the L-shape
(e.g. a node at `(-2.19, 0.22)`, and with some seeds astronomically far, up to
`~1e23`). Detected with the new `assess_nodes_outside_domain(Node, SDF)`.

**Root cause (two compounding bugs in the Python port):**

1. **`PolyMshr_CntrdPly` misaligned centroids with seeds.** It looped over *all*
   regions and `continue`d on any degenerate region (containing scipy's `-1`
   infinite-vertex marker or `<3` vertices) while only advancing a shared
   `counter`. When a seed cell was skipped, a later *reflected-point* region
   filled its slot, so `centroids[i]` no longer belonged to seed `i`. The seeds
   got scrambled, pushing one outside the domain -> its cell became unbounded ->
   more degenerate regions -> runaway divergence. (MATLAB loops `el = 1:NElem`
   one-to-one and never skips.)

2. **Corrupt cells launched seeds to infinity.** Even bounded cells can have
   far-away vertices when generators are near-collinear (aggravated by closely
   spaced fixed points), producing a centroid far outside the domain.

**Fix (`src/pyPolyMesher/pyPolyMesher.py`):**

- Rewrote `PolyMshr_CntrdPly(Element, Node, NElem, P=None)` to iterate exactly
  over the first `NElem` seed regions, one centroid per seed (aligned), dropping
  `-1` markers and using only finite vertices. Degenerate cells fall back to the
  seed's previous position `P[el]` (never to far vertices).
- Added a **stability guard** in the Lloyd loop: any centroid that is non-finite
  or lands outside the domain (`SDF > Tol`) is rejected and the seed is kept in
  place, so the iteration recovers.

**Result (12 random seeds, NElem=100):**

| Config | worst `max_outside` before | after |
|---|---|---|
| no PFix | 8.97e12 | 0.0588 |
| convex only | 2.9e23 | 0.1239 |
| apex + convex | 4.4e16 | 0.1230 |
| apex + convex + reentrant | 1.6e23 | 0.1204 |

Catastrophic divergence eliminated; all 62 tests still pass. Remaining outside
distances (~0.05-0.12) are small boundary overshoots (Voronoi vertex just past a
straight wall) - a separate, milder issue to address next.

---

## 8. Boundary overshoots (FIXED)

**Symptom:** after the divergence fix (§7), meshes still had a few nodes a small
distance *outside* straight walls (the free top/right edges of the L), `< 0.13`.
Validated over 10 random runs: avg 3.8 outside nodes/run, worst 0.125.
Cause: boundary Voronoi vertices shared by two boundary cells poke slightly past
a flat wall — inherent to discrete reflection.

**Fix (`src/pyPolyMesher/pyPolyMesher.py`):**
- Added `project_nodes_to_domain(Node, SDF, ...)`: pulls any node with `SDF > 0`
  back onto the boundary with Newton steps along the SDF gradient
  (`x <- x - d * grad(d)/|grad(d)|`). Interior nodes never move.
- `PolyMesher(..., snap_boundary=True)` now applies it as a final step (default on,
  exported for standalone use).

**Validation (10 randomized runs, now encoded in `tests/test_mesh_quality.py`):**

| | before snap | after snap |
|---|---|---|
| total outside nodes (10 runs) | 38 | **0** |
| worst `max_outside` | 0.1248 | **0.0000** |
| degenerate/sliver elements | 0 | **0** |

Generalizes to curved domains: Horn went 30 -> 0 outside nodes; Michell stayed 0.
Element quality preserved (min area > 0.014, no slivers). All 62 tests pass.
The 10-run grid was all clean.

### Always validate randomly
Because meshing is randomized, correctness must be checked over MANY runs, not
one. `tests/test_mesh_quality.py` runs the L-shape with multiple random seeds,
reports outside-node counts, and fails if any node falls outside the domain.

---

## 9. Boundary not exactly the original shape / persistent negative area error (FIXED)

**Symptom:** even with no nodes outside, the meshed geometry was slightly smaller
than the true domain (consistently **negative** area error, ~ -0.5%). The borders
were not exactly the original shape.

**Diagnosis (measured):**
- A pinned convex corner `(2,1)` had no exact node (nearest 0.0065 away) -> corner clipped.
- 7 of 47 boundary edges had their midpoint INSIDE the domain (up to 0.10 deep):
  a "reflection-gap notch" where a boundary cell's outer vertex is fixed by
  interior seeds, so the boundary dips inward and loses area.
- The earlier snap only pulled OUTSIDE nodes in; it never pushed these INSIDE
  boundary nodes back out to the wall, so the undershoot (negative area) remained.

**Fix (`src/pyPolyMesher/pyPolyMesher.py`):**
- `project_nodes_to_domain(..., node_ids=, max_step=)` generalized: given the
  boundary node ids it snaps each onto `d = 0` regardless of sign (the Newton
  step `x <- x - d*grad(d)/|grad(d)|` converges to the boundary from either
  side). `max_step` caps displacement (= element size) to prevent tangling.
- Added `_boundary_node_ids(Element, n)` (edges used by exactly one element).
- `PolyMesher(..., snap_boundary=True)` now snaps ALL boundary nodes onto the
  exact boundary, fixing both outside overshoots and inside notches.

**Validation (10 random runs, NElem=100):**

| | before | after |
|---|---|---|
| area error | ~ -0.485% (always negative) | mean |err| 0.027%, several runs exactly 0.000% |
| nodes outside | 0 | 0 |
| inverted/tangled elements | - | 0 |

Borders now follow the true L-shape outline exactly (verified in the randomized
`tests/test_mesh_quality.py` runs).
All 62 tests pass. Note: this boundary projection is a post-processing step added
on top of the MATLAB algorithm; the MATLAB reference does not project, so its
boundaries carry the same small reflection-gap deviations (just not measured).

---

## 10. Collinear boundary nodes introduced by snapping (FIXED)

**Symptom:** boundary elements sometimes had 3 nodes on a single straight wall
(two collinear edges within one element, e.g. a pentagon with 2 interior + 3
collinear boundary nodes).

**Cause:** the boundary projection in §9 snaps a reflection-gap notch vertex onto
the wall. That vertex then lies exactly on the straight line through its two
boundary neighbours (interior angle 180 deg), i.e. a redundant node creating two
collinear edges. Confirmed: snap=False -> 0 such nodes; snap=True (projection
only) -> 30 over 6 runs.

**Fix (`src/pyPolyMesher/pyPolyMesher.py`):**
- Added `PolyMshr_RmvCollinearNds(Node, Element)`: removes a node only if it is
  collinear (interior angle ~180 deg, and lying *between* its neighbours) in
  *every* element that contains it, then compacts the node list. The
  "collinear in all incident elements" rule keeps the mesh conforming (no
  hanging/T-junction nodes) and never touches true corners.
- Called right after the boundary projection in `PolyMesher`.

**Validation (10 random runs, correct per-element boundary detector):**

| | before | after |
|---|---|---|
| elements with 3 collinear boundary nodes | many (30 / 6 runs) | **0 / 10 runs** |
| area error | ~0 | ~0 (unchanged) |
| nodes outside / inverted / orphan | 0 / 0 / 0 | 0 / 0 / 0 |

(Note: a crude "any vertex near 180 deg" detector still reports ~2 hits/10 runs,
but those are *interior* 3-cell Voronoi junctions that are merely near-flat in
one of their cells -- a normal CVT feature, not a boundary defect.)
All 62 tests pass.
