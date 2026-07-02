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
