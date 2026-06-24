from pathlib import Path

import pytest

from pyPolyMesher.dxfImporter import dxf_polygon


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
