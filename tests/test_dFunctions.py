import numpy as np
import pytest

from pyPolyMesher.dFunctions import (
    dArc,
    dCircle,
    dDiff,
    dEllipse,
    dIntersect,
    dLine,
    dLineExact,
    dLoop,
    dPolygon,
    dRectangle,
    dRegularPolygon,
    dSpline,
    dStar,
    dUnion,
)


class TestSignedDistanceFields:
    def test_dLine(self):
        P = np.array([[0, 0], [1, 1], [2, 2], [0, 2]])
        result = dLine(P, 0, 0, 2, 2)

        assert result.shape == (4, 2)
        np.testing.assert_array_almost_equal(result[:3], np.zeros((3, 2)))
        np.testing.assert_array_almost_equal(
            result[3], np.array([-np.sqrt(2), -np.sqrt(2)])
        )

    def test_dLineExact(self):
        P = np.array([[0, 0], [1, 1], [2, 2], [0, 2], [2.1, 2]])
        result = dLineExact(P, 0, 0, 2, 2)

        assert result.shape == (5, 2)
        np.testing.assert_array_almost_equal(result[:3], np.zeros((3, 2)))
        np.testing.assert_array_almost_equal(
            result[3], np.array([-np.sqrt(2), -np.sqrt(2)])
        )
        np.testing.assert_array_almost_equal(result[4], np.array([0.1, 0.1]))

    def test_dCircle(self):
        P = np.array([[0, 0], [0, 1], [1, 1], [2, 2]])
        result = dCircle(P, 0, 0, 1)

        expected = np.array(
            [[-1, -1], [0, 0], [np.sqrt(2) - 1, np.sqrt(2) - 1], [np.sqrt(8) - 1, np.sqrt(8) - 1]]
        )
        np.testing.assert_array_almost_equal(result, expected)

    def test_dArc_full_circle_matches_dCircle(self):
        P = np.array([[0, 0], [0, 1], [1, 1], [2, 2]])
        arc = dArc(P, 0, 0, 1, 0.0, 2 * np.pi, ccw=True)
        circ = dCircle(P, 0, 0, 1)
        np.testing.assert_array_almost_equal(arc, circ)

    def test_dArc_quarter_on_arc_and_endpoints(self):
        # quarter circle from angle 0 to pi/2, center origin, r=2
        P = np.array([[2, 0], [0, 2], [2, 2], [1, 1], [3, 0]])
        result = dArc(P, 0, 0, 2, 0.0, np.pi / 2, ccw=True)
        assert result.shape == (5, 2)
        np.testing.assert_array_almost_equal(result[0, -1], 0.0)      # on arc
        np.testing.assert_array_almost_equal(result[1, -1], 0.0)      # on arc (endpoint)
        np.testing.assert_array_almost_equal(result[4, -1], 1.0)      # at start endpoint, r=2 -> dist 1 to (2,0)
        assert result[3, -1] < 0  # inside the full circle
        assert result[2, -1] > 0  # outside circle

    def test_dArc_outside_angular_span_uses_endpoints(self):
        # point at angle ~pi (opposite side) should be distance to nearest endpoint
        P = np.array([[-2, 0]])
        result = dArc(P, 0, 0, 2, 0.0, np.pi / 2, ccw=True)
        # nearest endpoint is (2,0) or (0,2); distance = sqrt(8)
        np.testing.assert_array_almost_equal(result[0, -1], np.sqrt(8))

    def test_dArc_clockwise_matches_counter_clockwise_magnitude(self):
        # cw arc from pi/2 to 0 traces the same quarter as ccw 0 to pi/2
        P = np.array([[2, 0], [0, 2], [2, 2], [1, 1], [3, 0], [-2, 0]])
        ccw = dArc(P, 0, 0, 2, 0.0, np.pi / 2, ccw=True)
        cw = dArc(P, 0, 0, 2, np.pi / 2, 0.0, ccw=False)
        np.testing.assert_array_almost_equal(np.abs(cw[:, -1]), np.abs(ccw[:, -1]))
        # signs must also match (same parent circle)
        np.testing.assert_array_almost_equal(cw, ccw)

    def test_dArc_center_point_and_offset_origin(self):
        # point at the arc center -> distance -r (inside parent circle)
        P = np.array([[5, 5]])
        result = dArc(P, 5, 5, 3, 0.0, np.pi / 2, ccw=True)
        np.testing.assert_array_almost_equal(result, np.array([[-3.0, -3.0]]))
        # offset start_angle and non-origin center: a point on the arc at angle pi/4
        P2 = np.array([[5 + 2 * np.cos(np.pi / 4), 5 + 2 * np.sin(np.pi / 4)]])
        r2 = dArc(P2, 5, 5, 2, np.pi / 4, np.pi, ccw=True)
        np.testing.assert_array_almost_equal(r2[:, -1], 0.0)

    def test_dSpline_open_degree1_is_straight_line(self):
        # degree-1 B-spline with 2 control points is a straight segment (0,0)-(2,0)
        P = np.array([[1, 0], [1, 1], [3, 0], [-1, 0]])
        result = dSpline(P, [(0, 0), (2, 0)], degree=1)
        assert result.shape == (4, 2)
        np.testing.assert_array_almost_equal(result[0, -1], 0.0)   # on the line
        np.testing.assert_array_almost_equal(result[1, -1], 1.0)   # 1 above
        np.testing.assert_array_almost_equal(result[2, -1], 1.0)   # 1 past end -> endpoint dist
        np.testing.assert_array_almost_equal(result[3, -1], 1.0)   # 1 before start

    def test_dSpline_open_clamped_passes_through_endpoints(self):
        # clamped cubic passes through first and last control point
        ctrl = [(0, 0), (1, 2), (2, 2), (3, 0)]
        P = np.array([[0, 0], [3, 0]])
        result = dSpline(P, ctrl, degree=3)
        np.testing.assert_array_almost_equal(result[:, -1], 0.0)

    def test_dSpline_closed_has_pip_sign(self):
        # closed quadratic spline approximating a loop around the origin
        ctrl = [(2, 0), (0, 2), (-2, 0), (0, -2)]
        P = np.array([[0, 0], [5, 5]])
        result = dSpline(P, ctrl, degree=2, closed=True)
        assert result[0, -1] < 0  # origin inside
        assert result[1, -1] > 0  # far outside

    def test_dLoop_four_lines_equals_dRectangle(self):
        segs = [
            {"type": "line", "p1": (0, 0), "p2": (2, 0)},
            {"type": "line", "p1": (2, 0), "p2": (2, 2)},
            {"type": "line", "p1": (2, 2), "p2": (0, 2)},
            {"type": "line", "p1": (0, 2), "p2": (0, 0)},
        ]
        # Exterior point is taken opposite an edge interior ([1,3], above the
        # top edge) rather than past a corner: dRectangle's max(x1-x, x-x2,
        # y1-y, y-y2) formula is exact for interior points and exterior points
        # that project onto an edge, but underestimates at exterior corners
        # (e.g. [3,3] -> 1 instead of the true sqrt(2)). dLoop computes the true
        # Euclidean per-segment distance, so the two agree only where dRectangle
        # is exact.
        P = np.array([[0, 0], [1, 1], [2, 2], [1, 3]])
        result = dLoop(P, segs)
        rect = dRectangle(P, 0, 2, 0, 2)
        assert result.shape == (4, 5)  # 4 segments + 1 signed
        np.testing.assert_array_almost_equal(result[:, -1], rect[:, -1])

    def test_dLoop_square_sign_inside_outside(self):
        segs = [
            {"type": "line", "p1": (1, 0), "p2": (1, 1)},
            {"type": "line", "p1": (1, 1), "p2": (0, 1)},
            {"type": "line", "p1": (0, 1), "p2": (0, 0)},
            {"type": "line", "p1": (0, 0), "p2": (1, 0)},
        ]
        P = np.array([[0.5, 0.5], [2, 2]])
        result = dLoop(P, segs)
        assert result[0, -1] < 0  # inside
        assert result[1, -1] > 0  # outside

    def test_dLoop_multi_loop_parity_hole(self):
        outer = [
            {"type": "line", "p1": (0, 0), "p2": (4, 0)},
            {"type": "line", "p1": (4, 0), "p2": (4, 4)},
            {"type": "line", "p1": (4, 4), "p2": (0, 4)},
            {"type": "line", "p1": (0, 4), "p2": (0, 0)},
        ]
        hole = [
            {"type": "line", "p1": (1, 1), "p2": (3, 1)},
            {"type": "line", "p1": (3, 1), "p2": (3, 3)},
            {"type": "line", "p1": (3, 3), "p2": (1, 3)},
            {"type": "line", "p1": (1, 3), "p2": (1, 1)},
        ]
        P = np.array([[2, 2], [0.5, 0.5], [5, 5]])
        result = dLoop(P, [outer, hole])
        assert result[0, -1] > 0   # in the hole -> outside domain (even parity)
        assert result[1, -1] < 0   # in the solid ring
        assert result[2, -1] > 0   # outside everything
        assert result.shape == (3, 9)  # 8 segments + 1 signed

    def test_dLoop_line_collinear_outside_uses_endpoint_distance(self):
        # point (3,0) is collinear with the bottom edge (0,0)-(2,0) but past its
        # endpoint; true distance is 1 (to corner (2,0)), NOT 0.
        segs = [
            {"type": "line", "p1": (0, 0), "p2": (2, 0)},
            {"type": "line", "p1": (2, 0), "p2": (2, 2)},
            {"type": "line", "p1": (2, 2), "p2": (0, 2)},
            {"type": "line", "p1": (0, 2), "p2": (0, 0)},
        ]
        P = np.array([[3, 0]])
        result = dLoop(P, segs)
        np.testing.assert_array_almost_equal(result[0, -1], 1.0)

    def test_dRectangle(self):
        P = np.array([[0, 0], [1, 1], [2, 2], [3, 3]])
        result = dRectangle(P, 0, 2, 0, 2)

        expected = np.array(
            [
                [0, -2, 0, -2, 0],
                [-1, -1, -1, -1, -1],
                [-2, 0, -2, 0, 0],
                [-3, 1, -3, 1, 1],
            ]
        )
        np.testing.assert_array_almost_equal(result, expected)

    def test_dPolygon(self):
        P = np.array([[0, 0], [1, 1], [2, 2], [3, 3]])
        vertices = [(0, 0), (2, 0), (2, 2), (0, 2)]
        result = dPolygon(P, vertices.copy())

        expected = np.array(
            [
                [0.0, -2.0, -2.0, 0.0, -0.0],
                [-1.0, -1.0, -1.0, -1.0, -1.0],
                [-2.0, 0.0, 0.0, -2.0, 0.0],
                [-3.162278, 1.414214, 1.414214, -3.162278, 1.414214],
            ]
        )
        np.testing.assert_array_almost_equal(result, expected)

    def test_dPolygon_reversed_vertex_order_and_closed_ring(self):
        P = np.array([[1, 1], [3, 3]])
        ccw_vertices = [(0, 0), (2, 0), (2, 2), (0, 2)]
        reversed_vertices = list(reversed(ccw_vertices))
        closed_vertices = ccw_vertices + [ccw_vertices[0]]

        ccw_result = dPolygon(P, ccw_vertices.copy())
        reversed_result = dPolygon(P, reversed_vertices.copy())
        closed_result = dPolygon(P, closed_vertices.copy())

        np.testing.assert_array_almost_equal(ccw_result[:, -1], reversed_result[:, -1])
        np.testing.assert_array_almost_equal(ccw_result[:, -1], closed_result[:, -1])
        assert ccw_result[0, -1] < 0  # inside
        assert ccw_result[1, -1] > 0  # outside

    def test_dEllipse(self):
        P = np.array([[0, 0], [2, 0], [0, 1], [3, 0]])
        result = dEllipse(P, 0, 0, 2, 1)

        assert result.shape == (4, 2)
        assert result[0, -1] < 0  # center is inside
        np.testing.assert_array_almost_equal(result[1], np.array([0.0, 0.0]))
        np.testing.assert_array_almost_equal(result[2], np.array([0.0, 0.0]))
        assert result[3, -1] > 0  # outside

    def test_dRegularPolygon(self):
        P = np.array([[0, 0], [1, 0], [2, 0]])
        result = dRegularPolygon(P, 0, 0, 1, 6)

        assert result.shape == (3, 7)
        assert result[0, -1] < 0  # center is inside
        assert abs(result[1, -1]) < 1e-12  # on a vertex for a 6-gon at theta=0
        assert result[2, -1] > 0  # outside

    def test_dStar(self):
        P = np.array([[0, 0], [1, 0], [2, 0]])
        result = dStar(P, 0, 0, 2, 1, 5)

        assert result.shape == (3, 11)
        assert result[0, -1] < 0  # center is inside the star
        assert result[1, -1] <= 0  # should be on or inside depending on orientation
        assert result[2, -1] >= 0  # outside or on boundary

    @pytest.fixture
    def point_grid(self):
        x = np.linspace(-5, 5, 11)
        y = np.linspace(-5, 5, 11)
        X, Y = np.meshgrid(x, y)
        return np.column_stack((X.ravel(), Y.ravel()))

    def test_dDiff(self, point_grid):
        circle = dCircle(point_grid, 0, 0, 3)
        rectangle = dRectangle(point_grid, -2, 2, -1, 1)

        result = dDiff(circle, rectangle)
        expected = np.maximum(circle[:, -1], -rectangle[:, -1])

        assert result.shape == (point_grid.shape[0], 6)
        np.testing.assert_array_almost_equal(result[:, -1], expected)
        assert result[np.where((point_grid[:, 0] == 0) & (point_grid[:, 1] == 0))[0][0], -1] > 0
        assert result[np.where((point_grid[:, 0] == 4) & (point_grid[:, 1] == 4))[0][0], -1] > 0

    def test_dIntersect(self, point_grid):
        circle1 = dCircle(point_grid, -1, 0, 2)
        circle2 = dCircle(point_grid, 1, 0, 2)

        result = dIntersect(circle1, circle2)
        expected = np.maximum(circle1[:, -1], circle2[:, -1])

        assert result.shape == (point_grid.shape[0], 3)
        np.testing.assert_array_almost_equal(result[:, -1], expected)
        assert result[np.where((point_grid[:, 0] == 0) & (point_grid[:, 1] == 0))[0][0], -1] < 0
        assert result[np.where((point_grid[:, 0] == -4) & (point_grid[:, 1] == 0))[0][0], -1] > 0

    def test_dUnion(self, point_grid):
        circle = dCircle(point_grid, 2, 2, 1.5)
        rectangle = dRectangle(point_grid, -3, 0, -1, 1)

        result = dUnion(circle, rectangle)
        expected = np.minimum(circle[:, -1], rectangle[:, -1])

        assert result.shape == (point_grid.shape[0], 6)
        np.testing.assert_array_almost_equal(result[:, -1], expected)
        assert result[np.where((point_grid[:, 0] == -2) & (point_grid[:, 1] == 0))[0][0], -1] < 0
        assert result[np.where((point_grid[:, 0] == 4) & (point_grid[:, 1] == 4))[0][0], -1] > 0

    def test_nested_boolean_composition(self, point_grid):
        outer = dCircle(point_grid, 0, 0, 4)
        hole = dCircle(point_grid, 0, 0, 1.5)
        bar = dRectangle(point_grid, -4, 4, -0.5, 0.5)

        intersected = dIntersect(outer, bar)
        nested = dDiff(intersected, hole)
        expected = np.maximum(intersected[:, -1], -hole[:, -1])

        assert nested.shape == (point_grid.shape[0], 7)
        np.testing.assert_array_almost_equal(nested[:, -1], expected)
        assert nested[np.where((point_grid[:, 0] == 0) & (point_grid[:, 1] == 0))[0][0], -1] > 0
        assert nested[np.where((point_grid[:, 0] == 0) & (point_grid[:, 1] == 4))[0][0], -1] > 0


if __name__ == "__main__":
    pytest.main()
