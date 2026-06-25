"""Tests for boundary tools: project_nodes_to_domain, _boundary_node_ids,
PolyMshr_RmvCollinearNds."""
import numpy as np
import pytest

from pyPolyMesher import project_nodes_to_domain
from pyPolyMesher.pyPolyMesher import _boundary_node_ids, PolyMshr_RmvCollinearNds
from pyPolyMesher.dFunctions import dRectangle, dCircle


def unit_square_sdf(P):
    return dRectangle(P, 0.0, 1.0, 0.0, 1.0)


def circle_sdf(P):
    return dCircle(P, 0.0, 0.0, 1.0)


# --------------------------- project_nodes_to_domain ---------------------------


def test_default_pulls_outside_nodes_onto_boundary():
    Node = np.array([[0.5, 0.5], [1.3, 0.5], [0.5, 1.4]])
    out = project_nodes_to_domain(Node, unit_square_sdf)
    d = unit_square_sdf(out)[:, -1]
    assert np.all(d <= 1e-8)                 # nothing left outside
    assert out[0] == pytest.approx([0.5, 0.5])  # interior node untouched
    assert out[1] == pytest.approx([1.0, 0.5])  # projected onto x = 1
    assert out[2] == pytest.approx([0.5, 1.0])  # projected onto y = 1


def test_node_ids_projects_inside_node_outward_to_boundary():
    # An inside boundary node should be pushed OUT onto the wall when targeted.
    Node = np.array([[0.5, 0.5], [0.9, 0.5]])  # second node 0.1 inside x = 1
    out = project_nodes_to_domain(Node, unit_square_sdf, node_ids=[1])
    assert out[1] == pytest.approx([1.0, 0.5], abs=1e-6)
    assert out[0] == pytest.approx([0.5, 0.5])  # untargeted node untouched


def test_projects_onto_curved_boundary():
    Node = np.array([[1.5, 0.0], [0.0, 0.7]])  # outside / inside the unit circle
    out = project_nodes_to_domain(Node, circle_sdf, node_ids=[0, 1])
    r = np.linalg.norm(out, axis=1)
    assert r == pytest.approx([1.0, 1.0], abs=1e-4)


def test_max_step_caps_displacement():
    # Node 1 at (0.8, 0.5) is 0.2 inside the x = 1 wall; projecting would move it
    # 0.2, but max_step caps the displacement at 0.1.
    Node = np.array([[0.5, 0.5], [0.8, 0.5]])
    out = project_nodes_to_domain(Node, unit_square_sdf, node_ids=[1], max_step=0.1)
    moved = np.linalg.norm(out[1] - Node[1])
    assert moved == pytest.approx(0.1, abs=1e-9)


# ----------------------------- _boundary_node_ids ------------------------------


def test_boundary_node_ids_two_squares_share_internal_edge():
    # Two unit squares sharing edge (1,2); the shared edge nodes are interior.
    Node = np.array(
        [
            [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0],  # 0..3 left square
            [2.0, 0.0], [2.0, 1.0],                          # 4,5 right square
        ]
    )
    Element = [np.array([0, 1, 2, 3]), np.array([1, 4, 5, 2])]
    bnd = set(_boundary_node_ids(Element, Node.shape[0]).tolist())
    # Every node lies on the outer boundary of the 2x1 rectangle here.
    assert bnd == {0, 1, 2, 3, 4, 5}


def test_boundary_node_ids_excludes_interior_node():
    # A central node shared by all elements is interior (no boundary edge once).
    # Build a square split into 4 triangles around center node 4.
    Node = np.array(
        [[0, 0], [1, 0], [1, 1], [0, 1], [0.5, 0.5]], dtype=float
    )
    Element = [
        np.array([0, 1, 4]),
        np.array([1, 2, 4]),
        np.array([2, 3, 4]),
        np.array([3, 0, 4]),
    ]
    bnd = set(_boundary_node_ids(Element, Node.shape[0]).tolist())
    assert 4 not in bnd
    assert bnd == {0, 1, 2, 3}


# --------------------------- PolyMshr_RmvCollinearNds --------------------------


def test_removes_collinear_node_shared_by_all_elements():
    # Node 2 (0.5, 0.0) lies on the straight edge between 1 and 3 in BOTH elems.
    Node = np.array(
        [
            [0.0, 0.0],   # 0
            [0.0, 0.0],   # placeholder (unused position kept simple) -> overwrite below
        ]
    )
    Node = np.array(
        [
            [0.0, 1.0],   # 0 interior-ish apex of elem A
            [0.0, 0.0],   # 1 corner
            [0.5, 0.0],   # 2 COLLINEAR on bottom edge
            [1.0, 0.0],   # 3 corner
            [1.0, 1.0],   # 4 apex of elem B
        ]
    )
    Element = [np.array([0, 1, 2, 3]), np.array([3, 2, 4])]  # both contain node 2
    # node 2 is collinear in elem A (1-2-3 straight). In elem B (3-2-4) it is the
    # vertex between (1,0)->(0.5,0) and (0.5,0)->(1,1): not collinear -> kept.
    NewNode, NewElement = PolyMshr_RmvCollinearNds(Node, Element)
    used = {int(i) for e in NewElement for i in e}
    # node 2 is NOT collinear in every element, so it must be preserved.
    assert NewNode.shape[0] == 5
    assert len(used) == 5


def test_removes_node_collinear_everywhere():
    # Node 1 lies on the straight line 0-2 and is the only "extra" point; both
    # elements see it as collinear, so it must be dropped and lists compacted.
    Node = np.array(
        [
            [0.0, 0.0],   # 0
            [0.5, 0.0],   # 1 collinear on bottom edge
            [1.0, 0.0],   # 2
            [0.5, 1.0],   # 3 apex
        ]
    )
    # Two triangles sharing the collinear chain on the bottom:
    Element = [np.array([0, 1, 3]), np.array([1, 2, 3])]
    # In elem A node 1 neighbors are 0 and 3 -> not collinear; so it stays.
    NewNode, NewElement = PolyMshr_RmvCollinearNds(Node, Element)
    assert NewNode.shape[0] == 4  # preserved (conforming rule)


def test_no_change_when_nothing_collinear():
    Node = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]])
    Element = [np.array([0, 1, 2])]
    NewNode, NewElement = PolyMshr_RmvCollinearNds(Node, Element)
    assert NewNode.shape[0] == 3
    assert [list(e) for e in NewElement] == [[0, 1, 2]]
