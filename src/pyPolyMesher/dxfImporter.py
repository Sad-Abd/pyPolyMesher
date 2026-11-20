from __future__ import annotations

from pathlib import Path
import warnings

import numpy as np

from pyPolyMesher.dFunctions import dLoop, _segment_points

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


def _polyline_vertices_to_segments(verts, closed=True):
    """Convert a list of (x, y, bulge) vertices into edge segments."""
    segs = []
    n = len(verts)
    for i in range(n):
        if i + 1 == n and not closed:
            break  # open polyline: no closing edge beyond last vertex
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
    if etype == "LWPOLYLINE":
        xs = [float(v) for v in data.get("10", [])]
        ys = [float(v) for v in data.get("20", [])]
        n = len(xs)
        bulges = [float(v) for v in data.get("42", [])]
        if len(bulges) < n:               # pad missing bulges with 0 (straight)
            bulges = bulges + [0.0] * (n - len(bulges))
        verts = list(zip(xs, ys, bulges))
        closed = bool(int(float(data.get("70", ["0"])[0])))  # bit 0 = closed flag
        return _polyline_vertices_to_segments(verts, closed=closed)
    if etype == "CIRCLE":
        return [{"type": "circle", "xc": _f(data, "10"), "yc": _f(data, "20"), "r": _f(data, "40")}]
    if etype == "ARC":
        return [{"type": "arc", "xc": _f(data, "10"), "yc": _f(data, "20"), "r": _f(data, "40"),
                 "start_angle": np.radians(_f(data, "50")), "end_angle": np.radians(_f(data, "51")), "ccw": True}]
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
        closed = (int(float(data.get("70", ["0"])[0])) & 1) != 0  # bit 0 = closed flag
        return [{"type": "spline", "ctrl": ctrl, "degree": degree, "knots": knots, "weights": weights, "closed": closed}]
    return []


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
    return a  # <0 for standard CCW, >0 for standard CW (this formula = -2 x standard signed area)


def _loop_polygon_dxf(loop):
    """Local alias for dFunctions._loop_polygon to avoid name shadowing."""
    from pyPolyMesher.dFunctions import _loop_polygon
    return _loop_polygon(loop)


def _loops_pfix(loops, bbox, angle_tol=15.0):
    # angle_tol must exceed the per-vertex turn of tessellated smooth curves
    # (a 64-sample circle/gon has turns ~360/63 ~ 5.7deg) so that smooth-curve
    # tessellation vertices are NOT mistaken for sharp corners. 15deg safely
    # filters those while still catching genuine sharp corners (>=15deg).
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
        ccw = _signed_area(poly) < 0  # standard convention: <0 = CCW
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
            closed = bool(int(float(data.get("70", ["0"])[0])))  # bit 0 = closed flag
            i += 1
            while i < len(entities) and entities[i][0] != "SEQEND":
                if entities[i][0] == "VERTEX":
                    vd = entities[i][1]
                    verts.append((_f(vd, "10"), _f(vd, "20"), _f(vd, "42", default="0")))
                i += 1
            i += 1  # skip SEQEND
            if verts:
                segments.extend(_polyline_vertices_to_segments(verts, closed=closed))
            continue
        segs = _entity_to_segments(etype, data)
        if segs:
            segments.extend(segs)
        i += 1

    if not segments:
        raise ValueError("No supported geometry entities found in DXF file.")
    return _chain_loops(segments)


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
