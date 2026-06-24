import numpy as np

from pyPolyMesher import Domain
from pyPolyMesher.dFunctions import dCircle, dRectangle


def test_domain_distfnc_delegates_to_configured_sdf():
    calls = []

    def sdf(points):
        calls.append(points.copy())
        return dRectangle(points, 0, 1, 0, 1)

    domain = Domain("unit-square", [0, 1, 0, 1], sdf)
    sample_points = np.array([[0.25, 0.25], [1.5, 1.5]])

    result = domain.DistFnc(sample_points)

    assert len(calls) == 1
    np.testing.assert_array_equal(calls[0], sample_points)
    np.testing.assert_array_equal(result, dRectangle(sample_points, 0, 1, 0, 1))


def test_domain_bndrycnds_returns_none_pair_when_no_callback_defined():
    domain = Domain("no-bc", [0, 1, 0, 1], lambda points: dRectangle(points, 0, 1, 0, 1))
    node = np.array([[0.0, 0.0], [1.0, 1.0]])

    result = domain.BndryCnds(node)

    assert result == [None, None]


def test_domain_bndrycnds_delegates_to_callback_with_nodes_and_bbox():
    received = {}

    def bc(nodes, bbox):
        received["nodes"] = nodes.copy()
        received["bbox"] = list(bbox)
        supp = np.array([[0, 1, 1]])
        load = np.array([[1, 0, -1]])
        return [supp, load]

    domain = Domain("with-bc", [-1, 2, -2, 3], lambda points: dRectangle(points, 0, 1, 0, 1), bc)
    nodes = np.array([[0.0, 0.0], [1.0, 1.0]])

    supp, load = domain.BndryCnds(nodes)

    np.testing.assert_array_equal(received["nodes"], nodes)
    assert received["bbox"] == [-1, 2, -2, 3]
    np.testing.assert_array_equal(supp, np.array([[0, 1, 1]]))
    np.testing.assert_array_equal(load, np.array([[1, 0, -1]]))


def test_domain_calculatearea_returns_plausible_circle_area_under_fixed_seed():
    domain = Domain("unit-circle", [-1, 1, -1, 1], lambda points: dCircle(points, 0, 0, 1))

    area = domain.CalculateArea(n=50_000)

    assert 2.9 < area < 3.35
    assert abs(area - np.pi) < 0.2
