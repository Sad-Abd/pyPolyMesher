"""Mesh-quality / performance tests.

These go beyond "runs without error": they mesh real domains and assert
geometric correctness of the result -- no nodes outside the domain, accurate
area, no inverted elements, captured corners, and no redundant collinear
boundary nodes. They also verify the divergence guard and the snap_boundary
post-processing actually improve the mesh.
"""
import numpy as np
import pytest

from pyPolyMesher import (
    Domain,
    PolyMesher,
    reentrant_corner_seeds,
    assess_nodes_outside_domain,
)
from pyPolyMesher.dFunctions import dRectangle, dDiff


def lshape_sdf(P):
    return dDiff(dRectangle(P, 0, 2, 0, 2), dRectangle(P, 1, 2, 1, 2))


L_AREA = 3.0


def lshape_domain():
    h = np.sqrt(L_AREA / 100)
    PFix = [[1, 1], [0, 0], [2, 0], [2, 1], [1, 2], [0, 2]]
    PFix += reentrant_corner_seeds([1, 1], [1, 0], [0, 1], 1.6 * h)
    return Domain("L", [0, 2, 0, 2], lshape_sdf, PFix=PFix)


def signed_areas(Node, Element):
    out = []
    for e in Element:
        e = [int(i) for i in e]
        v = Node[e]
        x, y = v[:, 0], v[:, 1]
        out.append(0.5 * np.sum(x * np.roll(y, -1) - y * np.roll(x, -1)))
    return np.array(out)


def per_element_boundary_collinear(Node, Element, SDF, ang_tol=3.0):
    """Count element-vertices that are (a) on the boundary and (b) collinear
    with their two neighbours within the same element -- the '3 nodes on one
    wall' defect."""
    n = 0
    for e in Element:
        e = [int(i) for i in e]
        v = Node[e]
        m = len(e)
        for k in range(m):
            a, b, c = v[(k - 1) % m], v[k], v[(k + 1) % m]
            u, w = a - b, c - b
            nu, nw = np.linalg.norm(u), np.linalg.norm(w)
            if nu < 1e-12 or nw < 1e-12:
                continue
            ang = np.degrees(np.arccos(np.clip(np.dot(u, w) / (nu * nw), -1, 1)))
            if ang > 180 - ang_tol and abs(SDF(b.reshape(1, 2))[0, -1]) < 1e-6:
                n += 1
    return n


SEEDS = [0, 1, 2, 3, 4]


@pytest.mark.parametrize("seed", SEEDS)
def test_lshape_has_no_nodes_outside_domain(seed):
    np.random.seed(seed)
    Node, Element, *_ = PolyMesher(lshape_domain(), 100, 80)
    res = assess_nodes_outside_domain(Node, lshape_sdf, verbose=False)
    assert res["num_outside"] == 0, f"outside nodes: {res['num_outside']}"
    assert res["max_outside"] == pytest.approx(0.0, abs=1e-6)


@pytest.mark.parametrize("seed", SEEDS)
def test_lshape_area_is_accurate(seed):
    np.random.seed(seed)
    Node, Element, *_ = PolyMesher(lshape_domain(), 100, 80)
    total = float(np.sum(np.abs(signed_areas(Node, Element))))
    rel_err = abs(total - L_AREA) / L_AREA
    assert rel_err < 5e-3, f"area error {100 * rel_err:.3f}% (area={total:.5f})"


@pytest.mark.parametrize("seed", SEEDS)
def test_lshape_no_inverted_or_far_nodes(seed):
    np.random.seed(seed)
    Node, Element, *_ = PolyMesher(lshape_domain(), 100, 80)
    A = signed_areas(Node, Element)
    assert np.all(A > 0), "found non-positive (inverted/degenerate) element area"
    # No node should have escaped far from the bounding box (divergence guard).
    assert np.all(np.abs(Node) < 10.0)


@pytest.mark.parametrize("seed", SEEDS)
def test_lshape_no_collinear_boundary_nodes(seed):
    np.random.seed(seed)
    Node, Element, *_ = PolyMesher(lshape_domain(), 100, 80)
    defects = per_element_boundary_collinear(Node, Element, lshape_sdf)
    assert defects == 0, f"{defects} collinear boundary nodes within elements"


@pytest.mark.parametrize("seed", SEEDS)
def test_lshape_captures_reentrant_apex_and_corners(seed):
    np.random.seed(seed)
    Node, Element, *_ = PolyMesher(lshape_domain(), 100, 80)
    for corner in [[1, 1], [0, 0], [2, 0], [1, 2], [0, 2]]:
        gap = np.min(np.linalg.norm(Node - corner, axis=1))
        assert gap < 0.05, f"corner {corner} not captured (gap {gap:.4f})"


def test_snap_boundary_removes_overshoots_vs_unsnapped():
    np.random.seed(0)
    dom = lshape_domain()
    Node_no, *_ = PolyMesher(dom, 100, 80, snap_boundary=False)
    np.random.seed(0)
    Node_yes, *_ = PolyMesher(lshape_domain(), 100, 80, snap_boundary=True)
    out_no = assess_nodes_outside_domain(Node_no, lshape_sdf, verbose=False)
    out_yes = assess_nodes_outside_domain(Node_yes, lshape_sdf, verbose=False)
    assert out_yes["num_outside"] == 0
    assert out_yes["max_outside"] <= out_no["max_outside"]


def test_divergence_guard_keeps_seeds_bounded_without_fixed_points():
    # Without fixed points the bare algorithm used to occasionally fling a seed
    # to infinity; the centroid-alignment fix + guard must keep everything finite.
    np.random.seed(5)
    Node, Element, _, _, P = PolyMesher(
        Domain("L-nopfix", [0, 2, 0, 2], lshape_sdf), 100, 80
    )
    assert np.all(np.isfinite(Node))
    assert np.all(np.abs(P) < 10.0)
