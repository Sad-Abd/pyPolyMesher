import numpy as np
import pytest

from pyPolyMesher.dFunctions import (
    dCircle,
    dDiff,
    dEllipse,
    dIntersect,
    dLine,
    dLineExact,
    dPolygon,
    dRectangle,
    dRegularPolygon,
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
