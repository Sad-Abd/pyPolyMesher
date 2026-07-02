from pathlib import Path

import math
import pytest

from pyPolyMesher.dxfImporter import (
    dxf_polygon,
    dxf_extract_geometry,
    _entity_to_segments,
    _bulge_to_arc,
    _polyline_vertices_to_segments,
)


EXAMPLES_DIR = Path(__file__).resolve().parents[1] / "examples"


def test_example_dxf_assets_exist_and_are_readable():
    for name in ["circle.dxf", "polygon1.dxf", "__plg6.dxf"]:
        path = EXAMPLES_DIR / name
        assert path.is_file()
        assert path.read_text().strip()


def write_simple_line_dxf(path: Path, connected: bool = True):
    if connected:
        content = "\n".join(
            [
                "0", "LINE", "10", "0", "20", "0", "11", "1", "21", "0",
                "0", "LINE", "10", "1", "20", "0", "11", "1", "21", "1",
                "0", "LINE", "10", "1", "20", "1", "11", "0", "21", "1",
                "0", "LINE", "10", "0", "20", "1", "11", "0", "21", "0",
                "0", "EOF",
            ]
        )
    else:
        content = "\n".join(
            [
                "0", "LINE", "10", "0", "20", "0", "11", "1", "21", "0",
                "0", "LINE", "10", "2", "20", "0", "11", "2", "21", "1",
                "0", "EOF",
            ]
        )
    path.write_text(content)


def write_simple_lwpolyline_dxf(path: Path):
    path.write_text(
        "\n".join(
            [
                "0", "LWPOLYLINE", "10", "0", "20", "0", "10", "1", "20", "0", "10", "1", "20", "1", "10", "0", "20", "1",
                "0", "EOF",
            ]
        )
    )


def test_dxf_polygon_parses_line_entities_into_ordered_vertices(tmp_path):
    path = tmp_path / "lines.dxf"
    write_simple_line_dxf(path, connected=True)

    vertices = dxf_polygon(str(path))

    assert vertices == [(1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]


def test_dxf_polygon_parses_lwpolyline_fixture(tmp_path):
    path = tmp_path / "lwpolyline.dxf"
    write_simple_lwpolyline_dxf(path)

    vertices = dxf_polygon(str(path))

    assert vertices == [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]


def test_dxf_polygon_raises_for_disconnected_line_entities(tmp_path):
    path = tmp_path / "disconnected_lines.dxf"
    write_simple_line_dxf(path, connected=False)

    with pytest.raises(ValueError, match="Lines are not connected"):
        dxf_polygon(str(path))


def test_dxf_polygon_raises_for_mixed_entities(tmp_path):
    path = tmp_path / "mixed_entities.dxf"
    path.write_text(
        "\n".join(
            [
                "0", "LINE", "10", "0", "20", "0", "11", "1", "21", "0",
                "0", "LWPOLYLINE", "10", "1", "20", "0", "10", "1", "20", "1", "0", "EOF",
            ]
        )
    )

    with pytest.raises(ValueError, match="not all lines or all polylines"):
        dxf_polygon(str(path))


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
