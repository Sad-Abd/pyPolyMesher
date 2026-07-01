# Diverse DXF Geometry Import → mesh-ready Domain — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend `dxfImporter.py` beyond polygons so diverse hand-authored DXF geometries (CIRCLE, ARC, ELLIPSE, SPLINE, polyline bulge arcs, holes) import directly into a mesh-ready `Domain`, and add reusable exact-curve SDF primitives (`dArc`, `dSpline`, `dLoop`) to `dFunctions.py`.

**Architecture:** Keep and extend the existing **hand-rolled DXF tokenizer** (NOT `ezdxf` — it is un-importable in this environment due to a corrupt system font crashing its module-level font scan; the hand-rolled parser already works on hand-authored DXFs, which is the stated scope). Add three SDF primitives to `dFunctions.py` usable standalone. In `dxfImporter.py`, parse entities into typed segment dicts, chain them into closed loops, build a combined SDF via `dLoop` with containment-parity hole detection, auto-compute the bounding box, and return a `Domain`. A `tessellate_curves` toggle falls back to pure-polyline SDFs.

**Tech Stack:** Python 3.10, numpy, scipy, matplotlib, pytest. No new dependencies. Tests run via `.venv/bin/python -m pytest`.

**Baseline:** `101 passed, 1 failed` — the single failure is a pre-existing test referencing a missing `examples/circle.dxf` (Task 0 fixes it).

---

## File Structure

| File | Responsibility |
|---|---|
| `src/pyPolyMesher/dFunctions.py` | Add `dArc`, `dSpline` (+ de Boor helpers), `dLoop` (+ segment helpers `_segment_distance`, `_segment_points`, `_loop_polygon`). |
| `src/pyPolyMesher/dxfImporter.py` | Refactor tokenizer into `_read_tokens`/`_parse_entities`; add `_entity_to_segments`, `_polyline_vertices_to_segments`, `_bulge_to_arc`, `_chain_loops`, `_loops_bbox`, `_loops_pfix`, `_tessellate_segment`, `dxf_extract_geometry`, `dxf_to_domain`. Keep `dxf_polygon` backward-compatible. |
| `src/pyPolyMesher/pyPolyMesher.py` | Chunked `CalculateArea` (retains only signed column). |
| `src/pyPolyMesher/__init__.py` | Export `dxf_to_domain`, `dArc`, `dSpline`, `dLoop`. |
| `tests/test_dFunctions.py` | Tests for `dArc`, `dSpline`, `dLoop`. |
| `tests/test_dxfImporter.py` | Tests for new extractors + `dxf_to_domain`; existing tests stay green. |
| `examples/circle.dxf` | New hand-authored CIRCLE fixture. |
| `README.md` | Rewrite DXF section; add INSERT/block expansion to Roadmap. |

**Segment dict contract** (used by `dLoop`, `dxf_to_domain`, and all helpers):
```python
{"type": "line",     "p1": (x1, y1), "p2": (x2, y2)}
{"type": "arc",      "xc":, "yc":, "r":, "start_angle":, "end_angle":, "ccw": bool}
{"type": "circle",   "xc":, "yc":, "r":}
{"type": "ellipse",  "xc":, "yc":, "a":, "b":, "theta":}              # full ellipse
{"type": "ellipse_arc", "xc":, "yc":, "a":, "b":, "theta":, "start":, "end":}  # partial -> tessellated
{"type": "spline",   "ctrl": [(x,y),...], "degree":, "knots":, "weights":, "closed": bool}
{"type": "polyline", "points": [(x,y),...]}                          # tessellated edge (open or closed)
```

**SDF return-shape contract** (critical — `PolyMshr_Rflct` reads `NBdrySegs = d.shape[1]-1` and uses leading columns as per-segment distances, abs'd): every combined domain SDF returns `(N, S+1)` — leading `S` columns = absolute distance to each segment, last column = parity-signed combined distance. `dLoop` enforces this.

---

### Task 0: Fix baseline — add missing `examples/circle.dxf` fixture

**Files:**
- Create: `examples/circle.dxf`
- Test: `tests/test_dxfImporter.py::test_example_dxf_assets_exist_and_are_readable`

- [ ] **Step 1: Create the fixture**

Write `examples/circle.dxf` (a hand-authored minimal DXF with one CIRCLE entity, center (50,50), radius 25):

```
0
SECTION
2
ENTITIES
0
CIRCLE
10
50.0
20
50.0
40
25.0
0
ENDSEC
0
EOF
```

- [ ] **Step 2: Verify the previously-failing test now passes**

Run: `.venv/bin/python -m pytest tests/test_dxfImporter.py::test_example_dxf_assets_exist_and_are_readable -v`
Expected: PASS.

- [ ] **Step 3: Confirm full baseline is green**

Run: `.venv/bin/python -m pytest tests/test_dxfImporter.py -q`
Expected: `5 passed`.

- [ ] **Step 4: Commit**

```bash
git add examples/circle.dxf
git commit -m "test: add missing examples/circle.dxf fixture (fixes pre-existing baseline failure)"
```

---

### Task 1: `dArc` SDF primitive

**Files:**
- Modify: `src/pyPolyMesher/dFunctions.py` (add `dArc` after `dCircle`)
- Test: `tests/test_dFunctions.py`

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_dFunctions.py` (inside `TestSignedDistanceFields`, after `test_dCircle`):

```python
    def test_dArc_full_circle_matches_dCircle(self):
        P = np.array([[0, 0], [0, 1], [1, 1], [2, 2]])
        arc = dArc(P, 0, 0, 1, 0.0, 2 * np.pi, ccw=True)
        circ = dCircle(P, 0, 0, 1)
        np.testing.assert_array_almost_equal(arc, circ)

    def test_dArc_quarter_on_arc_and_endpoints(self):
        # quarter circle from angle 0 to pi/2, center origin, r=2
        P = np.array([[2, 0], [0, 2], [2, 2], [1, 1], [3, 0]])
        result = dArc(P, 0, 0, 2, 0.0, np.pi / 2, ccw=True)
        assert result.shape == (5, 2)
        np.testing.assert_array_almost_equal(result[0, -1], 0.0)      # on arc
        np.testing.assert_array_almost_equal(result[1, -1], 0.0)      # on arc (endpoint)
        np.testing.assert_array_almost_equal(result[4, -1], 1.0)      # at start endpoint, r=2 -> dist 1 to (2,0)
        assert result[3, -1] < 0  # inside the full circle
        assert result[2, -1] > 0  # outside circle

    def test_dArc_outside_angular_span_uses_endpoints(self):
        # point at angle ~pi (opposite side) should be distance to nearest endpoint
        P = np.array([[-2, 0]])
        result = dArc(P, 0, 0, 2, 0.0, np.pi / 2, ccw=True)
        # nearest endpoint is (2,0) or (0,2); distance = sqrt(8)
        np.testing.assert_array_almost_equal(result[0, -1], np.sqrt(8))
```

Add `dArc` to the import list at the top of the file (alphabetical, before `dCircle`).

- [ ] **Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest tests/test_dFunctions.py -k dArc -v`
Expected: FAIL (`ImportError: cannot import name 'dArc'`).

- [ ] **Step 3: Implement `dArc`**

Add to `src/pyPolyMesher/dFunctions.py` immediately after `dCircle`:

```python
def dArc(P, xc, yc, r, start_angle, end_angle, ccw=True):
    """
    Calculate the signed distance from points P to a circular arc.

    The magnitude is the distance to the arc curve: when a point's polar
    angle (about the arc center) lies within the arc's angular span, the
    distance is the radial distance to the circle; otherwise it is the
    distance to the nearest arc endpoint. The sign is negative inside the
    parent circle and positive outside (so a full arc reduces to ``dCircle``).

    Parameters:
        P (numpy.ndarray): Array of 2D points, shape (N, 2).
        xc, yc (float): Center of the parent circle.
        r (float): Radius.
        start_angle, end_angle (float): Arc endpoints in radians.
        ccw (bool): If True the arc runs counter-clockwise from start to end.

    Returns:
        numpy.ndarray: (N, 2) array of signed distances (both columns equal).
    """
    dx = P[:, 0] - xc
    dy = P[:, 1] - yc
    radial = np.sqrt(dx ** 2 + dy ** 2)
    theta = np.mod(np.arctan2(dy, dx), 2 * np.pi)

    ext = abs(end_angle - start_angle)
    full = ext >= 2 * np.pi - 1e-9
    if full:
        within = np.ones(P.shape[0], dtype=bool)
    else:
        s = np.mod(start_angle, 2 * np.pi)
        if ccw:
            span = np.mod(end_angle - start_angle, 2 * np.pi)
            rel = np.mod(theta - s, 2 * np.pi)
        else:
            span = np.mod(start_angle - end_angle, 2 * np.pi)
            rel = np.mod(s - theta, 2 * np.pi)
        within = rel <= span

    circ_dist = radial - r
    p_start = np.array([xc + r * np.cos(start_angle), yc + r * np.sin(start_angle)])
    p_end = np.array([xc + r * np.cos(end_angle), yc + r * np.sin(end_angle)])
    d_start = np.sqrt((P[:, 0] - p_start[0]) ** 2 + (P[:, 1] - p_start[1]) ** 2)
    d_end = np.sqrt((P[:, 0] - p_end[0]) ** 2 + (P[:, 1] - p_end[1]) ** 2)

    arc_dist = np.where(within, np.abs(circ_dist), np.minimum(d_start, d_end))
    signed = np.where(radial < r, -arc_dist, arc_dist)
    return np.column_stack((signed, signed))
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `.venv/bin/python -m pytest tests/test_dFunctions.py -k dArc -v`
Expected: 3 passed.

- [ ] **Step 5: Commit**

```bash
git add src/pyPolyMesher/dFunctions.py tests/test_dFunctions.py
git commit -m "feat(dFunctions): add exact dArc signed-distance primitive"
```

---

### Task 2: `dSpline` SDF primitive (adaptive B-spline tessellation)

**Files:**
- Modify: `src/pyPolyMesher/dFunctions.py` (add de Boor helpers + `dSpline` after `dArc`)
- Test: `tests/test_dFunctions.py`

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_dFunctions.py`:

```python
    def test_dSpline_open_degree1_is_straight_line(self):
        # degree-1 B-spline with 2 control points is a straight segment (0,0)-(2,0)
        P = np.array([[1, 0], [1, 1], [3, 0], [-1, 0]])
        result = dSpline(P, [(0, 0), (2, 0)], degree=1)
        assert result.shape == (4, 2)
        np.testing.assert_array_almost_equal(result[0, -1], 0.0)   # on the line
        np.testing.assert_array_almost_equal(result[1, -1], 1.0)   # 1 above
        np.testing.assert_array_almost_equal(result[2, -1], 1.0)   # 1 past end -> endpoint dist
        np.testing.assert_array_almost_equal(result[3, -1], 1.0)   # 1 before start

    def test_dSpline_open_clamped_passes_through_endpoints(self):
        # clamped cubic passes through first and last control point
        ctrl = [(0, 0), (1, 2), (2, 2), (3, 0)]
        P = np.array([[0, 0], [3, 0]])
        result = dSpline(P, ctrl, degree=3)
        np.testing.assert_array_almost_equal(result[:, -1], 0.0)

    def test_dSpline_closed_has_pip_sign(self):
        # closed quadratic spline approximating a loop around the origin
        ctrl = [(2, 0), (0, 2), (-2, 0), (0, -2)]
        P = np.array([[0, 0], [5, 5]])
        result = dSpline(P, ctrl, degree=2, closed=True)
        assert result[0, -1] < 0  # origin inside
        assert result[1, -1] > 0  # far outside
```

Add `dSpline` to the import list in the test file.

- [ ] **Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest tests/test_dFunctions.py -k dSpline -v`
Expected: FAIL (`ImportError: cannot import name 'dSpline'`).

- [ ] **Step 3: Implement de Boor helpers and `dSpline`**

Add to `src/pyPolyMesher/dFunctions.py` after `dArc`:

```python
def _build_clamped_knots(n_ctrl, degree):
    """Clamped uniform open knot vector for n_ctrl control points."""
    if n_ctrl <= degree:
        return np.array([0.0] * (degree + 1) + [1.0] * (degree + 1))
    internal = n_ctrl - degree - 1
    if internal <= 0:
        knots = [0.0] * (degree + 1) + [1.0] * (degree + 1)
    else:
        knots = [0.0] * (degree + 1) + list(np.linspace(0, 1, internal + 2)[1:-1]) + [1.0] * (degree + 1)
    return np.array(knots)


def _bspline_basis_span(knots, degree, u):
    """Index i of the knot span [u_i, u_{i+1}) containing u (clamped)."""
    n = len(knots) - degree - 1
    if u >= knots[n]:
        return n - 1
    if u <= knots[degree]:
        return degree
    low, high = degree, n
    while high - low > 1:
        mid = (low + high) // 2
        if knots[mid] <= u:
            low = mid
        else:
            high = mid
    return low


def _bspline_point_clamped(u, knots, ctrl, degree, weights):
    """Evaluate a (rational) B-spline point at u via de Boor in homogeneous coords."""
    i = _bspline_basis_span(knots, degree, u)
    p = degree
    if weights is not None:
        d = np.array([[ctrl[i - p + j][0] * weights[i - p + j],
                       ctrl[i - p + j][1] * weights[i - p + j],
                       weights[i - p + j]] for j in range(p + 1)], dtype=float)
    else:
        d = np.array([[ctrl[i - p + j][0], ctrl[i - p + j][1], 1.0] for j in range(p + 1)], dtype=float)
    for r in range(1, p + 1):
        for j in range(p, r - 1, -1):
            den = knots[i + 1 + j - r] - knots[i - p + j]
            alpha = 0.0 if den == 0 else (u - knots[i - p + j]) / den
            d[j] = (1 - alpha) * d[j - 1] + alpha * d[j]
    return np.array([d[p][0] / d[p][2], d[p][1] / d[p][2]])


def _bspline_point_periodic(u, ctrl, degree):
    """Evaluate a uniform periodic (closed) B-spline at parameter u."""
    n = len(ctrl)
    p = degree
    i = int(np.floor(u))
    if i >= n:
        i = n - 1
    d = [np.asarray(ctrl[(i - p + j) % n], dtype=float) for j in range(p + 1)]
    for r in range(1, p + 1):
        for j in range(p, r - 1, -1):
            den = p + 1 - r  # uniform unit-spaced knots
            a = (u - (i - p + j)) / den
            d[j] = (1 - a) * d[j - 1] + a * d[j]
    return d[p]


def _spline_samples(ctrl, degree, knots, weights, n_samples, closed):
    """Tessellate a B-spline into n_samples points."""
    ctrl = [np.asarray(c, dtype=float) for c in ctrl]
    if n_samples is None:
        n_samples = max(50, 10 * len(ctrl))
    if closed:
        u_lo, u_hi = degree - 1, len(ctrl)
        us = np.linspace(u_lo, u_hi, n_samples, endpoint=False)
        return np.array([_bspline_point_periodic(u, ctrl, degree) for u in us])
    if knots is None:
        knots = _build_clamped_knots(len(ctrl), degree)
    u_lo, u_hi = knots[degree], knots[len(knots) - degree - 1]
    us = np.linspace(u_lo, u_hi, n_samples)
    return np.array([_bspline_point_clamped(u, knots, ctrl, degree, weights) for u in us])


def dSpline(P, control_points, degree=3, knots=None, weights=None, n_samples=None, closed=False):
    """
    Signed distance from points P to a (rational) B-spline curve.

    The curve is adaptively tessellated into a polyline; the distance is the
    polyline distance (an honest approximation; ``n_samples`` may be increased
    for accuracy). For a closed spline the sign is computed via point-in-polygon
    over the tessellation (negative inside). For an open spline both columns hold
    the unsigned distance (the enclosing loop supplies the sign).

    Parameters:
        P (numpy.ndarray): Array of 2D points, shape (N, 2).
        control_points (sequence): Control points (x, y).
        degree (int): Spline degree (default 3).
        knots (array-like, optional): Knot vector (open splines). Auto-built if None.
        weights (array-like, optional): Rational weights. Ignored for closed splines.
        n_samples (int, optional): Tessellation resolution. Auto-chosen if None.
        closed (bool): Whether the spline is a closed loop.

    Returns:
        numpy.ndarray: (N, 2) array of signed distances (both columns equal).
    """
    P = np.asarray(P, dtype=float)
    pts = _spline_samples(control_points, degree, knots, weights, n_samples, closed)
    if closed:
        verts = [tuple(p) for p in pts]
        signed = dPolygon(P, list(verts))[:, -1]
        return np.column_stack((signed, signed))
    dist = np.full(P.shape[0], np.inf)
    for k in range(len(pts) - 1):
        d = dLineExact(P, pts[k][0], pts[k][1], pts[k + 1][0], pts[k + 1][1])[:, -1]
        dist = np.minimum(dist, np.abs(d))
    return np.column_stack((dist, dist))
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `.venv/bin/python -m pytest tests/test_dFunctions.py -k dSpline -v`
Expected: 3 passed.

- [ ] **Step 5: Commit**

```bash
git add src/pyPolyMesher/dFunctions.py tests/test_dFunctions.py
git commit -m "feat(dFunctions): add dSpline SDF via adaptive B-spline tessellation"
```

---

### Task 3: Segment helpers + `dLoop` SDF primitive

**Files:**
- Modify: `src/pyPolyMesher/dFunctions.py` (add `_segment_distance`, `_segment_points`, `_loop_polygon`, `dLoop` after `dSpline`)
- Test: `tests/test_dFunctions.py`

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_dFunctions.py`:

```python
    def test_dLoop_four_lines_equals_dRectangle(self):
        segs = [
            {"type": "line", "p1": (0, 0), "p2": (2, 0)},
            {"type": "line", "p1": (2, 0), "p2": (2, 2)},
            {"type": "line", "p1": (2, 2), "p2": (0, 2)},
            {"type": "line", "p1": (0, 2), "p2": (0, 0)},
        ]
        P = np.array([[0, 0], [1, 1], [2, 2], [3, 3]])
        result = dLoop(P, segs)
        rect = dRectangle(P, 0, 2, 0, 2)
        assert result.shape == (4, 5)  # 4 segments + 1 signed
        np.testing.assert_array_almost_equal(result[:, -1], rect[:, -1])

    def test_dLoop_square_sign_inside_outside(self):
        segs = [
            {"type": "line", "p1": (1, 0), "p2": (1, 1)},
            {"type": "line", "p1": (1, 1), "p2": (0, 1)},
            {"type": "line", "p1": (0, 1), "p2": (0, 0)},
            {"type": "line", "p1": (0, 0), "p2": (1, 0)},
        ]
        P = np.array([[0.5, 0.5], [2, 2]])
        result = dLoop(P, segs)
        assert result[0, -1] < 0  # inside
        assert result[1, -1] > 0  # outside

    def test_dLoop_multi_loop_parity_hole(self):
        # outer square 0..4, inner square hole 1..3 -> point in hole is OUTSIDE domain
        outer = [
            {"type": "line", "p1": (0, 0), "p2": (4, 0)},
            {"type": "line", "p1": (4, 0), "p2": (4, 4)},
            {"type": "line", "p1": (4, 4), "p2": (0, 4)},
            {"type": "line", "p1": (0, 4), "p2": (0, 0)},
        ]
        hole = [
            {"type": "line", "p1": (1, 1), "p2": (3, 1)},
            {"type": "line", "p1": (3, 1), "p2": (3, 3)},
            {"type": "line", "p1": (3, 3), "p2": (1, 3)},
            {"type": "line", "p1": (1, 3), "p2": (1, 1)},
        ]
        P = np.array([[2, 2], [0.5, 0.5], [5, 5]])
        result = dLoop(P, [outer, hole])
        assert result[0, -1] > 0   # in the hole -> outside domain (even parity)
        assert result[1, -1] < 0   # in the solid ring
        assert result[2, -1] > 0   # outside everything
        assert result.shape == (3, 9)  # 8 segments + 1 signed
```

Add `dLoop` to the import list in the test file.

- [ ] **Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest tests/test_dFunctions.py -k dLoop -v`
Expected: FAIL (`ImportError: cannot import name 'dLoop'`).

- [ ] **Step 3: Implement segment helpers and `dLoop`**

Add to `src/pyPolyMesher/dFunctions.py` after `dSpline`:

```python
def _segment_distance(P, seg):
    """Signed distance (N,) to a single segment dict (sign per primitive's convention)."""
    t = seg["type"]
    if t == "line":
        return dLineExact(P, seg["p1"][0], seg["p1"][1], seg["p2"][0], seg["p2"][1])[:, -1]
    if t == "arc":
        return dArc(P, seg["xc"], seg["yc"], seg["r"], seg["start_angle"], seg["end_angle"], seg.get("ccw", True))[:, -1]
    if t == "circle":
        return dCircle(P, seg["xc"], seg["yc"], seg["r"])[:, -1]
    if t == "ellipse":
        return dEllipse(P, seg["xc"], seg["yc"], seg["a"], seg["b"], seg.get("theta", 0.0))[:, -1]
    if t == "ellipse_arc":
        pts = np.asarray(_segment_points(seg, 64), dtype=float)
        dist = np.full(P.shape[0], np.inf)
        for k in range(len(pts) - 1):
            d = dLineExact(P, pts[k][0], pts[k][1], pts[k + 1][0], pts[k + 1][1])[:, -1]
            dist = np.minimum(dist, np.abs(d))
        return dist
    if t == "spline":
        return dSpline(P, seg["ctrl"], seg.get("degree", 3), seg.get("knots"), seg.get("weights"), seg.get("n_samples"), seg.get("closed", False))[:, -1]
    if t == "polyline":
        pts = np.asarray(seg["points"], dtype=float)
        dist = np.full(P.shape[0], np.inf)
        for k in range(len(pts) - 1):
            d = dLineExact(P, pts[k][0], pts[k][1], pts[k + 1][0], pts[k + 1][1])[:, -1]
            dist = np.minimum(dist, np.abs(d))
        # close the ring only if the stored points do not already repeat the first
        # point (a zero-length closing edge would make dLineExact return NaN and
        # corrupt the whole column via np.minimum).
        if seg.get("closed") and len(pts) > 2 and np.linalg.norm(pts[-1] - pts[0]) > 1e-9:
            d = dLineExact(P, pts[-1][0], pts[-1][1], pts[0][0], pts[0][1])[:, -1]
            dist = np.minimum(dist, np.abs(d))
        return dist
    raise ValueError(f"Unknown segment type: {t}")


def _segment_points(seg, n=64):
    """Sample points on a segment (for bbox / polygon / tessellation)."""
    t = seg["type"]
    if t == "line":
        return [np.asarray(seg["p1"], dtype=float), np.asarray(seg["p2"], dtype=float)]
    if t == "arc":
        a0, a1 = seg["start_angle"], seg["end_angle"]
        if not seg.get("ccw", True):
            a0, a1 = a1, a0
        ang = np.linspace(a0, a1, n)
        return np.column_stack((seg["xc"] + seg["r"] * np.cos(ang), seg["yc"] + seg["r"] * np.sin(ang)))
    if t == "circle":
        ang = np.linspace(0, 2 * np.pi, n)
        return np.column_stack((seg["xc"] + seg["r"] * np.cos(ang), seg["yc"] + seg["r"] * np.sin(ang)))
    if t == "ellipse":
        ang = np.linspace(0, 2 * np.pi, n)
        a, b, th = seg["a"], seg["b"], seg.get("theta", 0.0)
        x = seg["xc"] + a * np.cos(ang) * np.cos(th) - b * np.sin(ang) * np.sin(th)
        y = seg["yc"] + a * np.cos(ang) * np.sin(th) + b * np.sin(ang) * np.cos(th)
        return np.column_stack((x, y))
    if t == "ellipse_arc":
        ang = np.linspace(seg["start"], seg["end"], n)
        a, b, th = seg["a"], seg["b"], seg.get("theta", 0.0)
        x = seg["xc"] + a * np.cos(ang) * np.cos(th) - b * np.sin(ang) * np.sin(th)
        y = seg["yc"] + a * np.cos(ang) * np.sin(th) + b * np.sin(ang) * np.cos(th)
        return np.column_stack((x, y))
    if t == "spline":
        return _spline_samples(seg["ctrl"], seg.get("degree", 3), seg.get("knots"), seg.get("weights"), seg.get("n_samples", n), seg.get("closed", False))
    if t == "polyline":
        return np.asarray(seg["points"], dtype=float)
    raise ValueError(f"Unknown segment type: {t}")


def _loop_polygon(loop):
    """Tessellated polygon vertices of a closed loop (dedup consecutive duplicates)."""
    pts = []
    for seg in loop:
        sp = np.asarray(_segment_points(seg), dtype=float)
        if len(pts) == 0:
            pts.extend([tuple(p) for p in sp])
        else:
            if np.linalg.norm(np.asarray(sp[0]) - np.asarray(pts[-1])) < 1e-9:
                pts.extend([tuple(p) for p in sp[1:]])
            else:
                pts.extend([tuple(p) for p in sp])
    if len(pts) > 1 and np.linalg.norm(np.asarray(pts[0]) - np.asarray(pts[-1])) < 1e-9:
        pts = pts[:-1]
    return [tuple(p) for p in pts]


def dLoop(P, loops):
    """
    Signed distance from points P to one or more closed loops of mixed segments.

    ``loops`` is a list of loops; each loop is a list of segment dicts. For
    convenience, a single loop (a list of segment dicts) may be passed directly.
    The magnitude is the minimum distance to any segment of any loop. The sign is
    determined by containment parity: a point is inside the domain if it lies
    inside an odd number of loops (so nested loops form holes, islands, etc.).

    Returns:
        numpy.ndarray: (N, S+1) where S is the total number of segments; the
        first S columns are absolute per-segment distances and the last column is
        the parity-signed combined distance.
    """
    P = np.asarray(P, dtype=float)
    if loops and isinstance(loops[0], dict):
        loops = [loops]
    seg_cols = []
    loop_polys = []
    for loop in loops:
        for seg in loop:
            seg_cols.append(np.abs(_segment_distance(P, seg)))
        loop_polys.append(_loop_polygon(loop))
    if seg_cols:
        seg_array = np.column_stack(seg_cols)
        magnitude = seg_array.min(axis=1)
    else:
        seg_array = np.empty((P.shape[0], 0))
        magnitude = np.full(P.shape[0], np.inf)
    inside_count = np.zeros(P.shape[0], dtype=int)
    for poly in loop_polys:
        inside_count += _is_point_in_polygon(P, poly).astype(int)
    inside = (inside_count % 2 == 1)
    signed = magnitude * np.where(inside, -1.0, 1.0)
    return np.column_stack([seg_array, signed])
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `.venv/bin/python -m pytest tests/test_dFunctions.py -k "dArc or dSpline or dLoop" -v`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/pyPolyMesher/dFunctions.py tests/test_dFunctions.py
git commit -m "feat(dFunctions): add dLoop with segment helpers and parity hole detection"
```

---

### Task 4: Refactor dxfImporter tokenizer (`_read_tokens`, `_parse_entities`); keep `dxf_polygon` green

**Files:**
- Modify: `src/pyPolyMesher/dxfImporter.py`
- Test: `tests/test_dxfImporter.py` (existing tests must stay green)

- [ ] **Step 1: Add a regression test for `dxf_polygon` on an LWPOLYLINE with bulge**

Add to `tests/test_dxfImporter.py`:

```python
def test_dxf_polygon_lwpolyline_with_bulge_is_backward_compatible(tmp_path):
    # dxf_polygon still returns plain vertices (bulge ignored) for an LWPOLYLINE
    path = tmp_path / "lw.dxf"
    path.write_text(
        "\n".join([
            "0", "LWPOLYLINE",
            "10", "0", "20", "0", "42", "0.5",
            "10", "1", "20", "0", "42", "0",
            "10", "1", "20", "1", "42", "0",
            "10", "0", "20", "1", "42", "0",
            "0", "EOF",
        ])
    )
    vertices = dxf_polygon(str(path))
    assert vertices == [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
```

- [ ] **Step 2: Run tests to verify the new test fails (old parser mispairs with interleaved 42)**

Run: `.venv/bin/python -m pytest tests/test_dxfImporter.py -k bulge -v`
Expected: FAIL (old `_pairwise(sum(zip(xs, ys), ()))` still works because `42` values go into `data["42"]`, not `xs`/`ys` — so this may already PASS. If it passes, that is fine; the refactor in Step 3 must keep it passing.)

- [ ] **Step 3: Refactor the parser**

Replace the **entire contents** of `src/pyPolyMesher/dxfImporter.py` with:

```python
from __future__ import annotations

from pathlib import Path

import numpy as np

from pyPolyMesher.dFunctions import dLoop

_STRUCTURAL = {
    "SECTION", "ENDSEC", "EOF", "TABLE", "ENDTAB", "BLOCK", "ENDBLK",
    "HEADER", "CLASSES", "OBJECTS", "FOODS", "THUMBNAILIMAGE", "ACDSDATA",
    "APPID", "BLOCK_RECORD", "DICTIONARY", "DIMSTYLE", "LAYER", "LTYPE", "STYLE", "VPORT",
}


def _read_tokens(file_path):
    """Read a DXF file and return its stripped lines."""
    return [line.strip() for line in Path(file_path).read_text().splitlines()]


def _parse_entities(tokens):
    """
    Parse a DXF token list into (entity_type, {code: [values]}) pairs.

    Every group-code-0 entity that is not a structural/table record is returned,
    including VERTEX/SEQEND records (so POLYLINE grouping can be done by callers).
    """
    entities = []
    i = 0
    current = None
    current_data = {}

    def flush():
        nonlocal current, current_data
        if current is not None:
            entities.append((current, current_data))
        current = None
        current_data = {}

    while i < len(tokens) - 1:
        code = tokens[i].strip()
        value = tokens[i + 1].strip()
        i += 2
        if code != "0":
            if current is None:
                continue
            current_data.setdefault(code, []).append(value)
            continue
        if value in _STRUCTURAL:
            flush()
            continue
        flush()
        current = value

    flush()
    return entities


def dxf_polygon(file_path):
    """
    Extract polygon vertices from LINE, LWPOLYLINE, or POLYLINE DXF entities.

    (Backward-compatible legacy API.) For richer geometry use ``dxf_to_domain``.

    Parameters:
        file_path (str): Path to the DXF file.

    Returns:
        list: List of (x, y) vertex tuples.

    Raises:
        ValueError: If lines are not connected, multiple polylines exist, or
            entities are not all lines or all polylines.
    """
    tokens = _read_tokens(file_path)
    entities = _parse_entities(tokens)
    ents_type = [e[0] for e in entities]

    def _pairwise(iterable):
        it = iter(iterable)
        return zip(it, it)

    if entities and all(e == "LINE" for e in ents_type):
        vertices = []
        line_vertices = []
        for _, data in entities:
            x1 = float(data.get("10", ["0"])[0])
            y1 = float(data.get("20", ["0"])[0])
            x2 = float(data.get("11", ["0"])[0])
            y2 = float(data.get("21", ["0"])[0])
            line_vertices.append(((x1, y1), (x2, y2)))
        for (p0_start, p0_end), (p1_start, _) in zip(line_vertices, line_vertices[1:] + line_vertices[:1]):
            if p0_end == p1_start:
                vertices.append(p0_end)
            else:
                raise ValueError("Lines are not connected!")
        return vertices

    if entities and all(e == "LWPOLYLINE" for e in ents_type):
        if len(entities) > 1:
            raise ValueError("Does not support multiple POLYLINE.")
        _, data = entities[0]
        xs = data.get("10", [])
        ys = data.get("20", [])
        return [(float(x), float(y)) for x, y in _pairwise(sum(zip(xs, ys), ()))]

    if entities and all(e == "POLYLINE" for e in ents_type):
        if len(entities) > 1:
            raise ValueError("Does not support multiple POLYLINE.")
        vertices = []
        inside_polyline = False
        current_vertex = None
        for idx in range(0, len(tokens) - 1, 2):
            code = tokens[idx].strip()
            value = tokens[idx + 1].strip()
            if code == "0" and value == "POLYLINE":
                inside_polyline = True
                continue
            if inside_polyline and code == "0" and value == "VERTEX":
                current_vertex = {}
                continue
            if inside_polyline and code == "0" and value == "SEQEND":
                break
            if inside_polyline and current_vertex is not None and code in {"10", "20"}:
                current_vertex.setdefault(code, []).append(value)
                if code == "20" and "10" in current_vertex:
                    vertices.append((float(current_vertex["10"][-1]), float(current_vertex["20"][-1])))
                    current_vertex = None
        return vertices

    raise ValueError("The entities in the file are not all lines or all polylines.")
```

- [ ] **Step 4: Run the full dxfImporter test suite**

Run: `.venv/bin/python -m pytest tests/test_dxfImporter.py -v`
Expected: all pass (existing 4 + new bulge-backward-compat).

- [ ] **Step 5: Commit**

```bash
git add src/pyPolyMesher/dxfImporter.py tests/test_dxfImporter.py
git commit -m "refactor(dxfImporter): extract _read_tokens/_parse_entities; keep dxf_polygon green"
```

---

### Task 5: Entity → segment conversion + bulge → arc

**Files:**
- Modify: `src/pyPolyMesher/dxfImporter.py` (add `_f`, `_bulge_to_arc`, `_polyline_vertices_to_segments`, `_entity_to_segments`)
- Test: `tests/test_dxfImporter.py`

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_dxfImporter.py`:

```python
import math
from pyPolyMesher.dxfImporter import _entity_to_segments, _bulge_to_arc, _polyline_vertices_to_segments


def test_bulge_to_arc_half_circle():
    # bulge = 1 -> included angle 180 deg; chord (0,0)->(2,0), bulge>0 -> ccw, center (1,0), r=1
    seg = _bulge_to_arc((0.0, 0.0), (2.0, 0.0), 1.0)
    assert seg["type"] == "arc"
    assert math.isclose(seg["xc"], 1.0, abs_tol=1e-9)
    assert math.isclose(seg["yc"], 0.0, abs_tol=1e-9)
    assert math.isclose(seg["r"], 1.0, abs_tol=1e-9)
    assert seg["ccw"] is True


def test_bulge_zero_is_line():
    seg = _bulge_to_arc((0.0, 0.0), (2.0, 0.0), 0.0)
    assert seg["type"] == "line"


def test_entity_to_segments_circle():
    segs = _entity_to_segments("CIRCLE", {"10": ["5"], "20": ["5"], "40": ["2"]})
    assert segs == [{"type": "circle", "xc": 5.0, "yc": 5.0, "r": 2.0}]


def test_entity_to_segments_arc():
    segs = _entity_to_segments("ARC", {"10": ["0"], "20": ["0"], "40": ["2"], "50": ["0"], "51": ["1.5707963267948966"]})
    assert segs[0]["type"] == "arc"
    assert math.isclose(segs[0]["start_angle"], 0.0)
    assert math.isclose(segs[0]["end_angle"], math.pi / 2)
    assert segs[0]["ccw"] is True


def test_polyline_vertices_to_segments_with_bulge():
    verts = [(0.0, 0.0, 0.0), (2.0, 0.0, 1.0), (2.0, 2.0, 0.0), (0.0, 2.0, 0.0)]
    segs = _polyline_vertices_to_segments(verts)
    assert segs[0]["type"] == "line"
    assert segs[1]["type"] == "arc"
    assert segs[2]["type"] == "line"
    assert segs[3]["type"] == "line"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest tests/test_dxfImporter.py -k "bulge or entity_to_segments or polyline_vertices" -v`
Expected: FAIL (`ImportError`).

- [ ] **Step 3: Implement the converters**

Add to `src/pyPolyMesher/dxfImporter.py` (after `_parse_entities`, before `dxf_polygon`):

```python
def _f(data, code, default="0", index=0):
    """Float of the index-th value of a group code from parsed entity data."""
    vals = data.get(code, [default])
    return float(vals[index] if index < len(vals) else vals[-1])


def _bulge_to_arc(p1, p2, bulge):
    """Convert a polyline edge with a bulge value into a line or arc segment."""
    if abs(bulge) < 1e-12:
        return {"type": "line", "p1": (float(p1[0]), float(p1[1])), "p2": (float(p2[0]), float(p2[1]))}
    p1 = np.asarray(p1, dtype=float)
    p2 = np.asarray(p2, dtype=float)
    chord = p2 - p1
    L = np.linalg.norm(chord)
    alpha = 4.0 * np.arctan(abs(bulge))      # included angle
    r = L / (2.0 * np.sin(alpha / 2.0))
    mid = (p1 + p2) / 2.0
    perp = np.array([-chord[1], chord[0]]) / L   # left of p1->p2
    s = 1.0 if bulge > 0 else -1.0
    d = r * np.cos(alpha / 2.0)
    center = mid + s * perp * d
    start_angle = np.arctan2(p1[1] - center[1], p1[0] - center[0])
    end_angle = np.arctan2(p2[1] - center[1], p2[0] - center[0])
    return {"type": "arc", "xc": float(center[0]), "yc": float(center[1]), "r": r,
            "start_angle": float(start_angle), "end_angle": float(end_angle), "ccw": bulge > 0}


def _polyline_vertices_to_segments(verts):
    """Convert a list of (x, y, bulge) vertices into edge segments (closed ring)."""
    segs = []
    n = len(verts)
    for i in range(n):
        x1, y1, b1 = verts[i]
        x2, y2, _ = verts[(i + 1) % n]
        segs.append(_bulge_to_arc((x1, y1), (x2, y2), b1))
    return segs


def _entity_to_segments(etype, data):
    """Convert a single parsed entity (non-POLYLINE) into a list of segment dicts."""
    if etype == "LINE":
        x1, y1 = _f(data, "10"), _f(data, "20")
        x2, y2 = _f(data, "11"), _f(data, "21")
        return [{"type": "line", "p1": (x1, y1), "p2": (x2, y2)}]
    if etype == "CIRCLE":
        return [{"type": "circle", "xc": _f(data, "10"), "yc": _f(data, "20"), "r": _f(data, "40")}]
    if etype == "ARC":
        return [{"type": "arc", "xc": _f(data, "10"), "yc": _f(data, "20"), "r": _f(data, "40"),
                 "start_angle": _f(data, "50"), "end_angle": _f(data, "51"), "ccw": True}]
    if etype == "ELLIPSE":
        xc, yc = _f(data, "10"), _f(data, "20")
        mx, my = _f(data, "11"), _f(data, "21")          # major-axis endpoint vector
        ratio = _f(data, "40")                            # minor/major ratio
        a = np.hypot(mx, my)
        b = a * ratio
        theta = np.arctan2(my, mx)
        start = _f(data, "41", default="0")
        end = _f(data, "42", default=str(2 * np.pi))
        if abs(end - start) >= 2 * np.pi - 1e-9:
            return [{"type": "ellipse", "xc": xc, "yc": yc, "a": a, "b": b, "theta": float(theta)}]
        return [{"type": "ellipse_arc", "xc": xc, "yc": yc, "a": a, "b": b, "theta": float(theta),
                 "start": float(start), "end": float(end)}]
    if etype == "SPLINE":
        degree = int(float(data.get("71", ["3"])[0]))
        xs = [float(v) for v in data.get("10", [])]
        ys = [float(v) for v in data.get("20", [])]
        ctrl = list(zip(xs, ys))
        knots = [float(v) for v in data.get("40", [])] or None
        weights = [float(v) for v in data.get("41", [])] or None
        closed = bool(int(float(data.get("73", ["0"])[0])))
        return [{"type": "spline", "ctrl": ctrl, "degree": degree, "knots": knots, "weights": weights, "closed": closed}]
    return []
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `.venv/bin/python -m pytest tests/test_dxfImporter.py -k "bulge or entity_to_segments or polyline_vertices" -v`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/pyPolyMesher/dxfImporter.py tests/test_dxfImporter.py
git commit -m "feat(dxfImporter): entity->segment conversion with bulge->arc"
```

---

### Task 6: Loop chaining (`_chain_loops`) + `dxf_extract_geometry`

**Files:**
- Modify: `src/pyPolyMesher/dxfImporter.py` (add `_segment_endpoints`, `_is_self_closed`, `_flip_segment`, `_chain_loops`, `dxf_extract_geometry`)
- Test: `tests/test_dxfImporter.py`

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_dxfImporter.py`:

```python
from pyPolyMesher.dxfImporter import dxf_extract_geometry


def test_dxf_extract_geometry_circle_is_single_self_closed_loop(tmp_path):
    path = tmp_path / "c.dxf"
    path.write_text("\n".join([
        "0", "SECTION", "2", "ENTITIES",
        "0", "CIRCLE", "10", "0", "20", "0", "40", "5",
        "0", "ENDSEC", "0", "EOF",
    ]))
    loops = dxf_extract_geometry(str(path))
    assert len(loops) == 1
    assert len(loops[0]) == 1
    assert loops[0][0]["type"] == "circle"


def test_dxf_extract_geometry_chains_lines_into_one_loop(tmp_path):
    path = tmp_path / "sq.dxf"
    path.write_text("\n".join([
        "0", "SECTION", "2", "ENTITIES",
        "0", "LINE", "10", "0", "20", "0", "11", "2", "21", "0",
        "0", "LINE", "10", "2", "20", "0", "11", "2", "21", "2",
        "0", "LINE", "10", "2", "20", "2", "11", "0", "21", "2",
        "0", "LINE", "10", "0", "20", "2", "11", "0", "21", "0",
        "0", "ENDSEC", "0", "EOF",
    ]))
    loops = dxf_extract_geometry(str(path))
    assert len(loops) == 1
    assert len(loops[0]) == 4


def test_dxf_extract_geometry_insert_warns(tmp_path):
    path = tmp_path / "ins.dxf"
    path.write_text("\n".join([
        "0", "SECTION", "2", "ENTITIES",
        "0", "INSERT", "2", "MYBLOCK", "10", "1", "20", "1",
        "0", "CIRCLE", "10", "0", "20", "0", "40", "1",
        "0", "ENDSEC", "0", "EOF",
    ]))
    with pytest.warns(UserWarning, match="INSERT"):
        loops = dxf_extract_geometry(str(path))
    assert loops[0][0]["type"] == "circle"


def test_dxf_extract_geometry_unclosed_loop_raises(tmp_path):
    path = tmp_path / "open.dxf"
    path.write_text("\n".join([
        "0", "SECTION", "2", "ENTITIES",
        "0", "LINE", "10", "0", "20", "0", "11", "1", "21", "0",
        "0", "LINE", "10", "9", "20", "9", "11", "9", "21", "9",
        "0", "ENDSEC", "0", "EOF",
    ]))
    with pytest.raises(ValueError, match="Unclosed loop"):
        dxf_extract_geometry(str(path))


def test_dxf_extract_geometry_no_geometry_raises(tmp_path):
    path = tmp_path / "empty.dxf"
    path.write_text("\n".join(["0", "SECTION", "2", "ENTITIES", "0", "ENDSEC", "0", "EOF"]))
    with pytest.raises(ValueError, match="No supported geometry"):
        dxf_extract_geometry(str(path))
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest tests/test_dxfImporter.py -k "extract_geometry" -v`
Expected: FAIL (`ImportError`).

- [ ] **Step 3: Implement chaining + extractor**

Add to `src/pyPolyMesher/dxfImporter.py` (after `_entity_to_segments`):

```python
import warnings

from pyPolyMesher.dFunctions import _segment_points  # reused for endpoints/tessellation


def _segment_endpoints(seg):
    """Return (start, end) points of an open segment, or None if self-closed."""
    t = seg["type"]
    if t == "line":
        return (np.asarray(seg["p1"], dtype=float), np.asarray(seg["p2"], dtype=float))
    if t == "arc":
        s = np.array([seg["xc"] + seg["r"] * np.cos(seg["start_angle"]),
                      seg["yc"] + seg["r"] * np.sin(seg["start_angle"])])
        e = np.array([seg["xc"] + seg["r"] * np.cos(seg["end_angle"]),
                      seg["yc"] + seg["r"] * np.sin(seg["end_angle"])])
        return (s, e)
    if t == "ellipse_arc":
        pts = np.asarray(_segment_points(seg, 2), dtype=float)
        return (pts[0], pts[-1])
    if t == "spline":
        if seg.get("closed"):
            return None
        pts = np.asarray(_segment_points(seg, 2), dtype=float)
        return (pts[0], pts[-1])
    if t == "polyline":
        pts = np.asarray(seg["points"], dtype=float)
        return (pts[0], pts[-1])
    return None  # circle, full ellipse -> self-closed


def _is_self_closed(seg, tol=1e-6):
    if seg["type"] in {"circle", "ellipse"}:
        return True
    if seg.get("closed"):
        return True
    ep = _segment_endpoints(seg)
    if ep is None:
        return True
    return bool(np.linalg.norm(ep[0] - ep[1]) < tol)


def _flip_segment(seg):
    """Reverse the direction of an open segment."""
    t = seg["type"]
    if t == "line":
        return {"type": "line", "p1": seg["p2"], "p2": seg["p1"]}
    if t == "arc":
        return {**seg, "start_angle": seg["end_angle"], "end_angle": seg["start_angle"], "ccw": not seg.get("ccw", True)}
    if t == "ellipse_arc":
        return {**seg, "start": seg["end"], "end": seg["start"]}
    if t == "polyline":
        return {**seg, "points": list(reversed(seg["points"]))}
    if t == "spline":
        return {**seg, "ctrl": list(reversed(seg["ctrl"]))}
    return seg


def _chain_loops(segments, tol=1e-6):
    """Group open segments into closed loops; self-closed segments become 1-segment loops."""
    loops = []
    open_segs = []
    for seg in segments:
        if _is_self_closed(seg, tol):
            loops.append([seg])
        else:
            open_segs.append(seg)

    used = [False] * len(open_segs)
    for start in range(len(open_segs)):
        if used[start]:
            continue
        used[start] = True
        loop = [open_segs[start]]
        _, cur_end = _segment_endpoints(open_segs[start])
        loop_start = _segment_endpoints(open_segs[start])[0]
        while True:
            found = None
            flip = False
            for j in range(len(open_segs)):
                if used[j]:
                    continue
                s, e = _segment_endpoints(open_segs[j])
                if np.linalg.norm(s - cur_end) < tol:
                    found, flip = j, False
                    break
                if np.linalg.norm(e - cur_end) < tol:
                    found, flip = j, True
                    break
            if found is None:
                raise ValueError(f"Unclosed loop: no segment connects to point {tuple(cur_end)}.")
            seg = open_segs[found]
            if flip:
                seg = _flip_segment(seg)
            used[found] = True
            loop.append(seg)
            _, cur_end = _segment_endpoints(seg)
            if np.linalg.norm(cur_end - loop_start) < tol:
                break
        loops.append(loop)
    return loops


_NON_GEOMETRY = {"HATCH", "TEXT", "MTEXT", "POINT", "DIMENSION", "LEADER", "SOLID",
                 "3DFACE", "ATTRIB", "ATTDEF", "VERTEX", "SEQEND"}


def dxf_extract_geometry(file_path):
    """
    Parse a DXF file and return geometry as a list of closed loops of segments.

    Parameters:
        file_path (str): Path to the DXF file.

    Returns:
        list: ``[loop, ...]`` where each loop is a list of segment dicts.

    Raises:
        ValueError: If no supported geometry is found or a loop cannot be closed.
    Warns:
        UserWarning: If INSERT (block reference) entities are present (skipped).
    """
    tokens = _read_tokens(file_path)
    entities = _parse_entities(tokens)
    segments = []
    i = 0
    while i < len(entities):
        etype, data = entities[i]
        if etype == "INSERT":
            warnings.warn(
                "INSERT (block reference) entities are not supported and will be "
                "skipped. See README roadmap (block expansion is future work).",
                UserWarning,
            )
            i += 1
            continue
        if etype in _NON_GEOMETRY:
            i += 1
            continue
        if etype == "POLYLINE":
            verts = []
            i += 1
            while i < len(entities) and entities[i][0] != "SEQEND":
                if entities[i][0] == "VERTEX":
                    vd = entities[i][1]
                    verts.append((_f(vd, "10"), _f(vd, "20"), _f(vd, "42", default="0")))
                i += 1
            i += 1  # skip SEQEND
            if verts:
                segments.extend(_polyline_vertices_to_segments(verts))
            continue
        segs = _entity_to_segments(etype, data)
        if segs:
            segments.extend(segs)
        i += 1

    if not segments:
        raise ValueError("No supported geometry entities found in DXF file.")
    return _chain_loops(segments)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `.venv/bin/python -m pytest tests/test_dxfImporter.py -v`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/pyPolyMesher/dxfImporter.py tests/test_dxfImporter.py
git commit -m "feat(dxfImporter): loop chaining + dxf_extract_geometry"
```

---

### Task 7: Bounding box + fixed-point helpers

**Files:**
- Modify: `src/pyPolyMesher/dxfImporter.py` (add `_loops_bbox`, `_signed_area`, `_loop_polygon_dxf`, `_loops_pfix`)
- Test: `tests/test_dxfImporter.py`

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_dxfImporter.py`:

```python
from pyPolyMesher.dxfImporter import _loops_bbox, _loops_pfix


def test_loops_bbox_unit_square():
    loop = [
        {"type": "line", "p1": (0, 0), "p2": (1, 0)},
        {"type": "line", "p1": (1, 0), "p2": (1, 1)},
        {"type": "line", "p1": (1, 1), "p2": (0, 1)},
        {"type": "line", "p1": (0, 1), "p2": (0, 0)},
    ]
    bbox = _loops_bbox([loop])
    assert bbox[0] <= 0 and bbox[1] >= 1
    assert bbox[2] <= 0 and bbox[3] >= 1


def test_loops_bbox_circle():
    loop = [{"type": "circle", "xc": 5, "yc": 5, "r": 2}]
    bbox = _loops_bbox([loop])
    assert bbox[0] <= 3 and bbox[1] >= 7
    assert bbox[2] <= 3 and bbox[3] >= 7


def test_loops_pfix_l_shape_has_reentrant_seeds():
    # L-shape polygon vertices (CW): (0,0),(0,3),(1,3),(1,1),(3,1),(3,0)
    pts = [(0, 0), (0, 3), (1, 3), (1, 1), (3, 1), (3, 0)]
    loop = [{"type": "polyline", "points": pts + [pts[0]], "closed": True}]
    pfix = _loops_pfix([loop], [0, 3, 0, 3])
    assert len(pfix) >= 3  # re-entrant corner at (1,1) contributes 3 seeds
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest tests/test_dxfImporter.py -k "loops_bbox or loops_pfix" -v`
Expected: FAIL (`ImportError`).

- [ ] **Step 3: Implement bbox + pfix helpers**

Add to `src/pyPolyMesher/dxfImporter.py` (after `_chain_loops`):

```python
def _loops_bbox(loops, pad_frac=0.02):
    """Auto bounding box [xmin, xmax, ymin, ymax] from all loop geometry, padded."""
    xs, ys = [], []
    for loop in loops:
        for seg in loop:
            pts = np.asarray(_segment_points(seg), dtype=float)
            xs.extend(pts[:, 0].tolist())
            ys.extend(pts[:, 1].tolist())
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    pad = pad_frac * max(xmax - xmin, ymax - ymin, 1e-9)
    return [xmin - pad, xmax + pad, ymin - pad, ymax + pad]


def _signed_area(poly):
    a = 0.0
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        a += (x2 - x1) * (y2 + y1)
    return a  # >0 means CW per the codebase's _is_counter_clockwise convention


def _loop_polygon_dxf(loop):
    """Local alias for dFunctions._loop_polygon to avoid name shadowing."""
    from pyPolyMesher.dFunctions import _loop_polygon
    return _loop_polygon(loop)


def _loops_pfix(loops, bbox, angle_tol=5.0):
    """
    Fixed-point seeds at sharp loop corners. Convex sharp corners contribute the
    corner point itself; re-entrant (concave) corners contribute the three
    reentrant_corner_seeds from Talischi et al. (2012).
    """
    from pyPolyMesher.pyPolyMesher import reentrant_corner_seeds

    diag = np.hypot(bbox[1] - bbox[0], bbox[3] - bbox[2])
    delta = 0.03 * diag
    angle_tol_rad = np.radians(angle_tol)
    pfix = []
    for loop in loops:
        poly = [tuple(p) for p in _loop_polygon_dxf(loop)]
        if len(poly) < 3:
            continue
        ccw = _signed_area(poly) < 0  # codebase convention: sum>0 is CCW
        n = len(poly)
        for i in range(n):
            prev = np.asarray(poly[(i - 1) % n], dtype=float)
            cur = np.asarray(poly[i], dtype=float)
            nxt = np.asarray(poly[(i + 1) % n], dtype=float)
            e_in = cur - prev
            e_out = nxt - cur
            lin = np.linalg.norm(e_in)
            lout = np.linalg.norm(e_out)
            if lin < 1e-12 or lout < 1e-12:
                continue
            cos_t = np.dot(e_in, e_out) / (lin * lout)
            turn = np.arccos(np.clip(cos_t, -1.0, 1.0))
            if turn < angle_tol_rad:
                continue
            cross = e_in[0] * e_out[1] - e_in[1] * e_out[0]
            concave = (cross < 0) if ccw else (cross > 0)
            if concave:
                e1 = -e_in / lin       # wall direction back along incoming edge
                e2 = e_out / lout      # wall direction along outgoing edge
                pfix.extend(reentrant_corner_seeds(cur, e1, e2, delta))
            else:
                pfix.append([float(cur[0]), float(cur[1])])
    if not pfix:
        return []
    arr = np.unique(np.asarray(pfix, dtype=float), axis=0)
    return arr.tolist()
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `.venv/bin/python -m pytest tests/test_dxfImporter.py -k "loops_bbox or loops_pfix" -v`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/pyPolyMesher/dxfImporter.py tests/test_dxfImporter.py
git commit -m "feat(dxfImporter): auto bbox + sharp/re-entrant corner fixed points"
```

---

### Task 8: `dxf_to_domain` — assemble mesh-ready Domain

**Files:**
- Modify: `src/pyPolyMesher/dxfImporter.py` (add `_tessellate_segment`, `dxf_to_domain`)
- Test: `tests/test_dxfImporter.py`

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_dxfImporter.py`:

```python
import numpy as np
from pyPolyMesher.dxfImporter import dxf_to_domain


def test_dxf_to_domain_circle_sdf_and_bbox(tmp_path):
    path = tmp_path / "c.dxf"
    path.write_text("\n".join([
        "0", "SECTION", "2", "ENTITIES",
        "0", "CIRCLE", "10", "5", "20", "5", "40", "2",
        "0", "ENDSEC", "0", "EOF",
    ]))
    dom = dxf_to_domain(str(path))
    assert dom.name == "c"
    assert dom.BdBox[0] <= 3 and dom.BdBox[1] >= 7
    P = np.array([[5, 5], [5, 7], [10, 10]])
    sdf = dom.SDF(P)
    assert sdf[0, -1] < 0  # center inside
    assert abs(sdf[1, -1]) < 1e-9  # on boundary
    assert sdf[2, -1] > 0  # outside


def test_dxf_to_domain_rectangle_with_circular_hole(tmp_path):
    # outer 0..6 square, circular hole at (3,3) r=1
    path = tmp_path / "rh.dxf"
    path.write_text("\n".join([
        "0", "SECTION", "2", "ENTITIES",
        "0", "LINE", "10", "0", "20", "0", "11", "6", "21", "0",
        "0", "LINE", "10", "6", "20", "0", "11", "6", "21", "6",
        "0", "LINE", "10", "6", "20", "6", "11", "0", "21", "6",
        "0", "LINE", "10", "0", "20", "6", "11", "0", "21", "0",
        "0", "CIRCLE", "10", "3", "20", "3", "40", "1",
        "0", "ENDSEC", "0", "EOF",
    ]))
    dom = dxf_to_domain(str(path), detect_fixed_points=False)
    P = np.array([[3, 3], [0.5, 0.5], [7, 7]])
    sdf = dom.SDF(P)
    assert sdf[0, -1] > 0   # in the hole -> outside domain
    assert sdf[1, -1] < 0   # in the solid ring
    assert sdf[2, -1] > 0   # outside everything


def test_dxf_to_domain_lwpolyline_with_bulge_fillet(tmp_path):
    # square with a 90deg fillet written as LWPOLYLINE with a bulge
    path = tmp_path / "f.dxf"
    path.write_text("\n".join([
        "0", "SECTION", "2", "ENTITIES",
        "0", "LWPOLYLINE",
        "10", "0", "20", "0", "42", "0",
        "10", "2", "20", "0", "42", "0.4142135623730951",  # 90deg fillet
        "10", "2", "20", "2", "42", "0",
        "10", "0", "20", "2", "42", "0",
        "0", "ENDSEC", "0", "EOF",
    ]))
    dom = dxf_to_domain(str(path), detect_fixed_points=False)
    P = np.array([[1, 1], [5, 5]])
    sdf = dom.SDF(P)
    assert sdf[0, -1] < 0
    assert sdf[1, -1] > 0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest tests/test_dxfImporter.py -k "dxf_to_domain" -v`
Expected: FAIL (`ImportError`).

- [ ] **Step 3: Implement `dxf_to_domain`**

Add to `src/pyPolyMesher/dxfImporter.py` (after `_loops_pfix`):

```python
def _tessellate_segment(seg, n=64):
    """Convert any segment into a polyline segment (tessellate_curves mode)."""
    if seg["type"] == "polyline":
        return seg
    pts = np.asarray(_segment_points(seg, n), dtype=float)
    closed = seg["type"] in {"circle", "ellipse"} or seg.get("closed", False)
    if closed and len(pts) > 1 and np.linalg.norm(pts[0] - pts[-1]) > 1e-9:
        pts = np.vstack([pts, pts[0]])
    return {"type": "polyline", "points": [tuple(p) for p in pts], "closed": closed}


def dxf_to_domain(file_path, tessellate_curves=False, detect_fixed_points=True, name=None):
    """
    Import a DXF file as a mesh-ready :class:`pyPolyMesher.Domain`.

    Parses LINE, LWPOLYLINE, POLYLINE (with bulge arcs), CIRCLE, ARC, ELLIPSE,
    and SPLINE entities, chains them into closed loops, detects holes by
    containment parity, auto-computes the bounding box, and builds the signed
    distance function. INSERT (block reference) entities are skipped with a
    warning (hand-authored DXFs only; block expansion is future work).

    Parameters:
        file_path (str): Path to the DXF file.
        tessellate_curves (bool): If True, sample every curve into polylines and
            build the SDF from polygonal loops only (the robust fallback).
            Default False uses exact arc/circle/ellipse primitives.
        detect_fixed_points (bool): If True, seed ``Domain.PFix`` at sharp and
            re-entrant corners to improve mesh quality. Default True.
        name (str, optional): Domain name. Defaults to the file stem.

    Returns:
        pyPolyMesher.Domain: A domain ready for ``PolyMesher(domain, ...)``.

    Raises:
        ValueError: If no supported geometry is found or a loop cannot be closed.
    Warns:
        UserWarning: If INSERT entities are present.
    """
    from pyPolyMesher.pyPolyMesher import Domain

    loops = dxf_extract_geometry(file_path)
    if tessellate_curves:
        loops = [[_tessellate_segment(seg) for seg in loop] for loop in loops]
    bbox = _loops_bbox(loops)

    def SDF(P):
        return dLoop(P, loops)

    pfix = _loops_pfix(loops, bbox) if detect_fixed_points else []
    dom_name = name if name is not None else Path(file_path).stem
    return Domain(dom_name, bbox, SDF, PFix=pfix)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `.venv/bin/python -m pytest tests/test_dxfImporter.py -v`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/pyPolyMesher/dxfImporter.py tests/test_dxfImporter.py
git commit -m "feat(dxfImporter): dxf_to_domain assembles mesh-ready Domain"
```

---

### Task 9: `tessellate_curves` mode equivalence test

**Files:**
- Test: `tests/test_dxfImporter.py`

- [ ] **Step 1: Write the equivalence test**

Add to `tests/test_dxfImporter.py`:

```python
def test_dxf_to_domain_tessellate_matches_exact_within_tol(tmp_path):
    path = tmp_path / "c.dxf"
    path.write_text("\n".join([
        "0", "SECTION", "2", "ENTITIES",
        "0", "CIRCLE", "10", "5", "20", "5", "40", "2",
        "0", "ENDSEC", "0", "EOF",
    ]))
    dom_exact = dxf_to_domain(str(path), tessellate_curves=False, detect_fixed_points=False)
    dom_tess = dxf_to_domain(str(path), tessellate_curves=True, detect_fixed_points=False)
    P = np.array([[5, 5], [5, 7], [6, 6], [10, 10]])
    e = dom_exact.SDF(P)[:, -1]
    t = dom_tess.SDF(P)[:, -1]
    assert np.all(np.sign(e) == np.sign(t))
    np.testing.assert_allclose(np.abs(t), np.abs(e), atol=0.2)
```

- [ ] **Step 2: Run the test**

Run: `.venv/bin/python -m pytest tests/test_dxfImporter.py -k tessellate_matches -v`
Expected: PASS.

If it fails on magnitude tolerance, bump `_tessellate_segment`'s `n` default to 128 and re-run.

- [ ] **Step 3: Commit**

```bash
git add tests/test_dxfImporter.py
git commit -m "test(dxfImporter): tessellate_curves mode matches exact within tolerance"
```

---

### Task 10: Chunked `CalculateArea`

**Files:**
- Modify: `src/pyPolyMesher/pyPolyMesher.py` (`Domain.CalculateArea`)
- Test: `tests/test_domain.py`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_domain.py` (ensure `import pytest` and `pytest.approx` are available):

```python
    def test_calculate_area_chunked_for_multi_segment_domain(self):
        from pyPolyMesher.exampleDomains import MbbDomain
        # Mbb is a 3x1 rectangle; chunked path must still give ~3.0
        area = MbbDomain.CalculateArea(n=20000)
        assert area == pytest.approx(3.0, rel=0.05)
```

- [ ] **Step 2: Run test to verify current behavior**

Run: `.venv/bin/python -m pytest tests/test_domain.py -k calculate_area -v`
Expected: may PASS already; the chunked rewrite must preserve the value.

- [ ] **Step 3: Rewrite `CalculateArea` to be chunked**

In `src/pyPolyMesher/pyPolyMesher.py`, replace the body of `CalculateArea` (after its docstring) with:

```python
    def CalculateArea(self, n=1_000_000, chunk=100_000):
        """
        Calculates the approximate area of the domain using the Monte Carlo method,
        evaluating the SDF in chunks so that only the signed column is retained.
        This keeps memory bounded for domains with many boundary segments.

        Args:
            n (int, optional): Total random samples. Default 1,000,000.
            chunk (int, optional): Samples per SDF evaluation. Default 100,000.

        Returns:
            float: The approximate area of the domain.
        """
        xmin, xmax, ymin, ymax = self.BdBox
        total_area = (xmax - xmin) * (ymax - ymin)
        points_inside = 0
        remaining = n
        while remaining > 0:
            k = min(chunk, remaining)
            pts = np.random.uniform(low=[xmin, ymin], high=[xmax, ymax], size=(k, 2))
            sdf_values = self.SDF(pts)
            points_inside += int(np.sum(sdf_values[:, -1] <= 0))
            remaining -= k
        return (points_inside / n) * total_area
```

- [ ] **Step 4: Run the domain tests**

Run: `.venv/bin/python -m pytest tests/test_domain.py -v`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/pyPolyMesher/pyPolyMesher.py tests/test_domain.py
git commit -m "perf(Domain): chunked CalculateArea retains only signed column"
```

---

### Task 11: Public exports

**Files:**
- Modify: `src/pyPolyMesher/__init__.py`
- Test: `tests/test_dxfImporter.py`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_dxfImporter.py`:

```python
def test_public_exports_available():
    import pyPolyMesher as ppm
    from pyPolyMesher import dxf_to_domain, dArc, dSpline, dLoop
    assert callable(ppm.dxf_to_domain)
    assert callable(ppm.dArc)
    assert callable(ppm.dSpline)
    assert callable(ppm.dLoop)
```

- [ ] **Step 2: Run to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_dxfImporter.py -k public_exports -v`
Expected: FAIL (`ImportError`).

- [ ] **Step 3: Update exports**

Replace `src/pyPolyMesher/__init__.py` with:

```python
from .pyPolyMesher import (
    PolyMesher,
    Domain,
    mesh_assessment,
    assess_nodes_outside_domain,
    reentrant_corner_seeds,
    project_nodes_to_domain,
)
from .dFunctions import dArc, dSpline, dLoop
from .dxfImporter import dxf_to_domain
```

- [ ] **Step 4: Run to verify it passes**

Run: `.venv/bin/python -m pytest tests/test_dxfImporter.py -k public_exports -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/pyPolyMesher/__init__.py tests/test_dxfImporter.py
git commit -m "feat: export dxf_to_domain, dArc, dSpline, dLoop"
```

---

### Task 12: PolyMesher smoke test on a DXF domain

**Files:**
- Test: `tests/test_dxfImporter.py`

- [ ] **Step 1: Write the smoke test**

Add to `tests/test_dxfImporter.py`:

```python
def test_dxf_to_domain_meshes_without_nodes_outside(tmp_path):
    from pyPolyMesher import PolyMesher, assess_nodes_outside_domain

    path = tmp_path / "rh.dxf"
    path.write_text("\n".join([
        "0", "SECTION", "2", "ENTITIES",
        "0", "LINE", "10", "0", "20", "0", "11", "6", "21", "0",
        "0", "LINE", "10", "6", "20", "0", "11", "6", "21", "6",
        "0", "LINE", "10", "6", "20", "6", "11", "0", "21", "6",
        "0", "LINE", "10", "0", "20", "6", "11", "0", "21", "0",
        "0", "CIRCLE", "10", "3", "20", "3", "40", "1",
        "0", "ENDSEC", "0", "EOF",
    ]))
    dom = dxf_to_domain(str(path), detect_fixed_points=False)
    Node, Element, Supp, Load, P = PolyMesher(dom, 30, 20, anim=False)
    report = assess_nodes_outside_domain(Node, dom.SDF, verbose=False)
    assert report["n_outside"] <= 1
```

- [ ] **Step 2: Run the smoke test**

Run: `.venv/bin/python -m pytest tests/test_dxfImporter.py -k meshes_without_nodes_outside -v`
Expected: PASS.

If it fails with many nodes outside, inspect `report` — likely the `dLoop` parity returns the wrong sign somewhere; debug the rectangle+hole case. Do not loosen the assertion to mask a real sign bug.

- [ ] **Step 3: Commit**

```bash
git add tests/test_dxfImporter.py
git commit -m "test(dxfImporter): PolyMesher smoke test on rectangle-with-hole DXF"
```

---

### Task 13: README + Roadmap update

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Update the DXF section**

In `README.md`, replace the "### Import Polygon Domain from DXF" section (heading + both code blocks + the `CalculateArea`/`mesh_assessment` paragraph) with:

````markdown
### Import Diverse Geometries from DXF

`dxf_to_domain` imports a hand-authored DXF (LINE, LWPOLYLINE/POLYLINE with bulge
arcs, CIRCLE, ARC, ELLIPSE, SPLINE) directly into a mesh-ready `Domain`. It
auto-computes the bounding box, detects holes by containment parity, and builds
the signed distance function from exact arc/circle/ellipse primitives (splines and
partial ellipses are tessellated). `INSERT` (block reference) entities are skipped
with a warning — block expansion is on the roadmap.

```python
from pyPolyMesher import PolyMesher, mesh_assessment, dxf_to_domain

dom = dxf_to_domain('examples/polygon1.dxf')
dom.Plot()

Node, Element, Supp, Load, P = PolyMesher(dom, 50, 100)
area = dom.CalculateArea()
metrics = mesh_assessment(Node, Element, area, verbose=True)
```

For a robust pure-polyline fallback (sample every curve), pass
`tessellate_curves=True`. The new SDF primitives `dArc`, `dSpline`, and `dLoop`
are also exported for building domains by hand without DXF.
````

Also update the table-of-contents link from `#import-polygon-domain-from-dxf` to `#import-diverse-geometries-from-dxf`.

- [ ] **Step 2: Add INSERT/block expansion to the Roadmap**

Find the `## Roadmap` section in `README.md` and add this item to its list:

```markdown
- [ ] DXF `INSERT`/block expansion: resolve block references with scale/rotation/translation (currently skipped with a warning)
```

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "docs: document dxf_to_domain; add block-expansion roadmap item"
```

---

### Task 14: Final full-suite verification

- [ ] **Step 1: Run the entire test suite**

Run: `.venv/bin/python -m pytest -q`
Expected: all pass (was `101 passed, 1 failed`; should now be ~120+ passed, 0 failed).

- [ ] **Step 2: If any test fails, debug with the systematic-debugging skill before proceeding**

Do not mark the plan complete with failing tests. Fix root causes, not symptoms.

- [ ] **Step 3: Confirm the example assets test and the new smoke tests all pass**

Run: `.venv/bin/python -m pytest tests/test_dxfImporter.py tests/test_dFunctions.py tests/test_domain.py -q`
Expected: all pass.
