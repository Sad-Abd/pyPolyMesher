from pathlib import Path
import runpy

import matplotlib.pyplot as plt
import pyPolyMesher.dxfImporter as dxf_importer


def test_examples_script_runs_without_gui_and_without_long_interaction(monkeypatch):
    repo_root = Path(__file__).resolve().parents[1]
    script = repo_root / "examples" / "Examples.py"
    polygon_dxf = repo_root / "examples" / "polygon1.dxf"

    assert script.is_file()
    assert polygon_dxf.is_file()
    assert polygon_dxf.read_text().strip()

    monkeypatch.setattr(plt, "show", lambda *args, **kwargs: None)
    monkeypatch.setattr(dxf_importer, "dxf_polygon", lambda path: [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
    monkeypatch.chdir(repo_root)

    namespace = runpy.run_path(str(script), run_name="__main__")

    assert "dxfDomain" in namespace
    assert namespace["dxfDomain"].name == "DXF Polygon Domain"
