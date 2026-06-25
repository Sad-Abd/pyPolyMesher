from __future__ import annotations

from pathlib import Path


def _pairwise(iterable):
    it = iter(iterable)
    return zip(it, it)


def dxf_polygon(file_path):
    """
    Extracts vertices from DXF entities of LINE, LWPOLYLINE, or POLYLINE types in the specified DXF file.

    Parameters:
        file_path (str): The path to the DXF file.

    Returns:
        list: A list of tuples representing the extracted vertices.

    Raises:
        ValueError: If the lines are not connected.
        ValueError: If there are multiple POLYLINE entities and the function does not support them.
        ValueError: If the entities in the file are neither all lines nor all polylines.
    """
    path = Path(file_path)
    tokens = path.read_text().splitlines()

    entities = []
    i = 0
    current = None
    current_data = {}

    def flush_current():
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
        # Group code 0 signals entity start; skip all structural markers.
        if value in {
            "SECTION", "ENDSEC", "EOF", "TABLE", "ENDTAB", "BLOCK", "ENDBLK",
            "HEADER", "CLASSES", "OBJECTS", "FOODS", "THUMBNAILIMAGE", "ACDSDATA",
        }:
            flush_current()
            continue
        # Skip non-geometry entities that appear in TABLES section.
        if value in {"APPID", "BLOCK_RECORD", "DICTIONARY", "DIMSTYLE", "LAYER", "LTYPE", "STYLE", "VPORT"}:
            flush_current()
            continue
        flush_current()
        current = value

    flush_current()

    ents = entities
    ents_type = [e[0] for e in ents]

    if ents and all(e == "LINE" for e in ents_type):
        vertices = []
        line_vertices = []
        for _, data in ents:
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

    elif ents and all(e == "LWPOLYLINE" for e in ents_type):
        if len(ents) > 1:
            raise ValueError("Does not support multiple POLYLINE.")
        _, data = ents[0]
        xs = data.get("10", [])
        ys = data.get("20", [])
        vertices = [(float(x), float(y)) for x, y in _pairwise(sum(zip(xs, ys), ()))]

    elif ents and all(e == "POLYLINE" for e in ents_type):
        if len(ents) > 1:
            raise ValueError("Does not support multiple POLYLINE.")
        # Parse VERTEX records following the POLYLINE entity.
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

    else:
        raise ValueError("The entities in the file are not all lines or all polylines.")

    return vertices
