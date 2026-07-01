# Diverse DXF Geometry Import → mesh-ready Domain

**Date:** 2026-07-02
**Status:** Approved
**Scope:** Extend `dxfImporter.py` beyond polygons so that diverse DXF geometries (circles, arcs, ellipses, splines, bulge arcs) can be imported and meshed directly. Add reusable exact-curve SDF primitives to `dFunctions.py`.

---

## Goal

Currently `dxfImporter.dxf_polygon()` only extracts polygon vertices from `LINE`/`LWPOLYLINE`/`POLYLINE` via a hand-rolled DXF tokenizer, and the user must manually wrap the result into a `Domain` with a hardcoded bounding box. We want:

1. Parse a broad set of DXF curve entities (circle, arc, ellipse, spline, polyline bulge arcs).
2. Produce a fully mesh-ready `Domain` (auto bounding box + signed distance function + hole detection) that can be passed straight to `PolyMesher(...)`.
3. Introduce first-class `dArc` / `dSpline` SDF primitives that are useful on their own — independent of DXF — for users building SDFs by hand.

## Decisions (from brainstorming)

| # | Question | Decision |
|---|----------|----------|
| Q1 | Importer output | `dxf_to_domain(path) -> Domain`, fully mesh-ready |
| Q2 | Curved entities | Hybrid exact curves + first-class `dArc`/`dSpline` primitives; user toggle to fall back to "tessellate everything to polylines" |
| Q3 | Hole classification | Pure geometric containment (containment-parity), no CAD-direction or layer-name assumptions |
| Q4 | `dSpline` fidelity | Exact `dArc`/`dCircle`/`dEllipse`; `dSpline` = adaptively-tessellated polyline distance behind a clean `dSpline(P, ...)` API (no NURBS closest-point solver) |
| Q5 | `INSERT` (blocks) | Hand-authored DXFs only — warn and skip `INSERT`; full block expansion is future work (README roadmap) |

## Architecture & Files

### `dFunctions.py` — new public SDF primitives

- **`dArc(P, xc, yc, r, start_angle, end_angle, ccw=True)`**
  Exact unsigned distance to a circular arc: project query points onto the parent circle, clamp the angle to the arc's angular span, handle endpoints (nearest endpoint when the projection falls outside the span). Sign relative to the parent circle (negative inside the full circle). For a full arc (span = 2π) it reduces to `dCircle`.

- **`dSpline(P, control_points, degree=3, knots=None, weights=None, n_samples=None, closed=False)`**
  Clean public API mirroring DXF `SPLINE` parameters. Internally: adaptive B-spline tessellation into a dense polyline, then polyline distance (reuse `dLineExact`/`dPolygon` machinery). Closed splines get point-in-polygon-based sign; open splines return unsigned distance in both columns (the enclosing loop supplies the sign). Honest about the approximation; `n_samples` can be user-tuned or auto-selected from control-point span.

- **`dLoop(P, segments)`**
  Reusable building block: signed distance to a closed loop of mixed segment kinds (`line`/`arc`/`ellipse`/`spline`/`circle`). Magnitude = min over all segments of each segment's exact unsigned distance. Sign = point-in-polygon over the loop's tessellated vertices (uniformly negative inside, positive outside). `segments` is a list of dicts: `{"type": "line"|"arc"|"circle"|"ellipse"|"spline", ...params, "tessellation": [...]}`. This is the core reusable closed-boundary primitive.

### `dxfImporter.py` — ezdxf-based extractors

- **`dxf_extract_geometry(path)`** — parse the DXF with `ezdxf`, walk `msp.modelspace()` entities, emit per-loop segment lists. Kept public for users who want raw geometry.

- **`dxf_to_domain(path, tessellate_curves=False, detect_fixed_points=True, ...) -> Domain`** — builds the `Domain`:
  - chains segments into closed loops within an endpoint tolerance (raises `ValueError` with diagnostics on unclosed loops),
  - computes the auto bounding box,
  - assembles the SDF (see below),
  - optionally seeds `PFix` at sharp/re-entrant corners (reuses `reentrant_corner_seeds` for concave corners).

- **`dxf_polygon(file_path)`** — kept for backward compatibility; delegates to the new extractor for ezdxf-valid files, with the legacy tokenizer retained as a fallback so existing bare-DXF tests stay green.

### `__init__.py`
Export `dxf_to_domain`, `dArc`, `dSpline`, `dLoop`.

### `README.md`
- Rewrite the DXF section to use `dxf_to_domain` (auto bbox, holes, curves) with an example covering a circle + hole + arc + spline.
- Add "Block / INSERT expansion" to the Roadmap as future work.

### `tests/`
- Extend `tests/test_dFunctions.py` with `dArc`, `dSpline`, `dLoop` cases.
- Extend `tests/test_dxfImporter.py` with `dxf_extract_geometry` / `dxf_to_domain` cases and an ezdxf-generated fixture in `examples/`.

## Entity Coverage

| DXF entity | Handling |
|---|---|
| `LINE` | exact via `dLineExact` |
| `CIRCLE` | exact via `dCircle` (its own closed loop) |
| `ARC` | exact via `dArc` |
| `ELLIPSE` (full) | exact via `dEllipse` |
| `ELLIPSE` (partial / arc) | tessellated edge (exact ellipse-arc SDF is out of scope) |
| `LWPOLYLINE` / `POLYLINE` | per edge; **bulge** values → `dArc` arc segments, otherwise `dLineExact` |
| `SPLINE` | `dSpline` (tessellated) |
| `INSERT` | **warn and skip** (future: block expansion) |
| `HATCH`/`TEXT`/`POINT`/etc. | skip silently; summary warning for unsupported curve types |

## SDF Assembly & Holes

### Loop chaining
Connect line/arc/spline endpoints within a tolerance into closed loops; raise `ValueError` listing any loop that cannot close.

### Hole detection — containment parity per query point (no pre-classification)
The combined domain SDF is computed directly from all loops:

- `magnitude(P)` = min over **all** loops of (min over that loop's segments of unsigned distance)
- `inside(P)` = **odd** number of loops whose tessellated polygon contains `P` (ray-casting containment)
- `signed(P)` = `magnitude(P) * (inside ? -1 : +1)`

This is purely geometric, handles disjoint outer regions, holes, and islands-in-holes-in-islands at any nesting depth, and requires no `dDiff` bookkeeping.

### SDF return shape — satisfies the reflection contract
`PolyMshr_Rflct` reads `NBdrySegs = d.shape[1] - 1` and uses the **leading columns as per-boundary-segment distances** (`np.abs(d[:, 0:NBdrySegs]) < Alpha`) to reflect points across the specific nearby segment. `dPolygon` satisfies this by returning one column per edge.

The new combined domain SDF therefore returns `(N, S+1)` where `S` = total boundary segments across all loops; leading columns = each segment's unsigned distance, last column = signed distance. A "segment" = one DXF edge (one LINE; each LWPOLYLINE edge including bulge arcs; each ARC/CIRCLE/ELLIPSE/SPLINE). This mirrors `dPolygon`'s per-edge columns so `PolyMshr_Rflct` works unchanged.

## Bounding Box & Fixed Points

- `BdBox` auto-computed from all geometry extents (arcs via angular span, ellipses via parametric extent, splines via samples), padded slightly. No more hardcoded `[0,100,0,100]`.
- `detect_fixed_points=True` (default on): seed `PFix` at sharp corners of the loops (reusing `reentrant_corner_seeds` for concave ones) to improve mesh quality on DXF geometries with sharp features.

## `tessellate_curves` toggle (Q2)

`dxf_to_domain(path, tessellate_curves=False)`:
- `False` (default) → exact `dArc`/`dCircle`/`dEllipse` + tessellated `dSpline`/partial-ellipse.
- `True` → sample every curve into polylines; SDF built from polygonal loops only (the "Option A" fallback). Same parity-based hole handling.

## Targeted Improvement: chunked `CalculateArea`

`Domain.CalculateArea` evaluates the full `(N, S+1)` SDF array for 1,000,000 points before slicing `[:, -1]`. For complex DXFs (many segments) this is memory-heavy. Add a **chunked** evaluation path in `CalculateArea` that only retains the signed column. Low-risk; benefits all multi-segment domains (including `dPolygon` ones).

## Error Handling

- Unreadable / missing file → clear error.
- No geometry → `ValueError`.
- Unclosed loop → `ValueError` with diagnostics.
- `ezdxf` parse failure → wrapped helpful message.
- `INSERT` present → `warnings.warn` (visible, non-fatal).

## Testing

- **`dArc`**: exact checks at on-arc, center, endpoints, outside angular span; equals `dCircle` for a full circle; FD gradient sanity.
- **`dSpline`**: degenerate-to-line case; distance at known curve points; closed-spline sign via PIP.
- **`dLoop`**: 4-line square equals `dRectangle`; loop with one arc edge vs tessellated `dPolygon` magnitude; sign parity.
- **`dxf_to_domain`**: ezdxf-written fixture (rectangle + circular hole + arc fillet + spline edge) → SDF sign correct at probe points (solid negative, hole positive, outside positive); bbox correct; `INSERT` warning via `pytest.warns`; `tessellate_curves=True` matches exact within tolerance; small `PolyMesher` smoke test (few elements / iterations) runs and nodes stay inside.

## Out of Scope / Future Work

- Full `INSERT`/block expansion with transforms (scale/rotate/translate, recursive) — noted in README Roadmap.
- Exact SDF for partial ellipse arcs.
- True NURBS closest-point solver for `dSpline` (currently tessellated).
