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
