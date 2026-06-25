"""Tests for node-vs-SDF assessment: assess_nodes_outside_domain."""
import numpy as np
import pytest

from pyPolyMesher import assess_nodes_outside_domain
from pyPolyMesher.dFunctions import dRectangle


def unit_square_sdf(P):
    return dRectangle(P, 0.0, 1.0, 0.0, 1.0)


def test_all_inside_reports_zero():
    Node = np.array([[0.25, 0.25], [0.5, 0.5], [0.75, 0.4]])
    res = assess_nodes_outside_domain(Node, unit_square_sdf, verbose=False)
    assert res["num_outside"] == 0
    assert res["frac_outside"] == 0.0
    assert res["max_outside"] == 0.0
    assert res["mean_outside"] == 0.0
    assert res["sum_outside"] == 0.0
    assert res["outside_ids"].size == 0
    assert res["distances"].shape == (3,)


def test_detects_outside_nodes_and_distances():
    # Two nodes outside: one 0.5 to the right, one 0.2 above.
    Node = np.array(
        [
            [0.5, 0.5],   # inside
            [1.5, 0.5],   # outside by 0.5 (x > 1)
            [0.5, 1.2],   # outside by 0.2 (y > 1)
        ]
    )
    res = assess_nodes_outside_domain(Node, unit_square_sdf, verbose=False)
    assert res["num_outside"] == 2
    assert set(res["outside_ids"].tolist()) == {1, 2}
    assert res["max_outside"] == pytest.approx(0.5)
    assert res["mean_outside"] == pytest.approx((0.5 + 0.2) / 2)
    assert res["sum_outside"] == pytest.approx(0.7)
    assert res["frac_outside"] == pytest.approx(2 / 3)


def test_tolerance_ignores_boundary_nodes():
    # Node exactly on the boundary (d == 0) must not count as outside.
    Node = np.array([[1.0, 0.5], [0.5, 1.0]])
    res = assess_nodes_outside_domain(Node, unit_square_sdf, verbose=False)
    assert res["num_outside"] == 0


def test_verbose_prints_summary(capsys):
    Node = np.array([[0.5, 0.5], [2.0, 0.5]])
    assess_nodes_outside_domain(Node, unit_square_sdf, verbose=True)
    out = capsys.readouterr().out
    assert "Nodes outside domain" in out
    assert "max  outside distance" in out
