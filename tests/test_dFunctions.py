import numpy as np
import pytest
from pyPolyMesher.dFunctions import (
    dLine,
    dLineExact,
    dCircle,
    dRectangle,
    dPolygon,
    dDiff,
    dIntersect,
    dUnion,
)


class TestSignedDistanceFields:

    def test_dLine(self):
        P = np.array([[0, 0], [1, 1], [2, 2], [0, 2]])
        x1, y1, x2, y2 = 0, 0, 2, 2
        expected = np.array([[0, 0], [0, 0], [0, 0], [-1.41421356, -1.41421356]])
        result = dLine(P, x1, y1, x2, y2)
        np.testing.assert_array_almost_equal(result, expected)

    def test_dLineExact(self):
        P = np.array([[0, 0], [1, 1], [2, 2], [0, 2], [2.1, 2]])
        x1, y1, x2, y2 = 0, 0, 2, 2
        expected = np.array(
            [[0, 0], [0, 0], [0, 0], [-1.41421356, -1.41421356], [0.1, 0.1]]
        )
        result = dLineExact(P, x1, y1, x2, y2)
        np.testing.assert_array_almost_equal(result, expected)

    def test_dCircle(self):
        P = np.array([[0, 0], [0, 1], [1, 1], [2, 2]])
        xc, yc, r = 0, 0, 1
        expected = np.array(
            [[-1, -1], [0, 0], [0.41421356, 0.41421356], [1.82842712, 1.82842712]]
        )
        result = dCircle(P, xc, yc, r)
        np.testing.assert_array_almost_equal(result, expected)

    def test_dRectangle(self):
        P = np.array([[0, 0], [1, 1], [2, 2], [3, 3]])
        x1, x2, y1, y2 = 0, 2, 0, 2
        expected = np.array(
            [
                [0, -2, 0, -2, 0],
                [-1, -1, -1, -1, -1],
                [-2, 0, -2, 0, 0],
                [-3, 1, -3, 1, 1],
            ]
        )
        result = dRectangle(P, x1, x2, y1, y2)
        np.testing.assert_array_almost_equal(result, expected)

    def test_dPolygon(self):
        P = np.array([[0, 0], [1, 1], [2, 2], [3, 3]])
        vertices = [(0, 0), (2, 0), (2, 2), (0, 2)]
        expected = np.array(
            [
                [0.0, -2.0, -2.0, 0.0, -0.0],
                [-1.0, -1.0, -1.0, -1.0, -1.0],
                [-2.0, 0.0, 0.0, -2.0, 0.0],
                [-3.162278, 1.414214, 1.414214, -3.162278, 1.414214],
            ]
        )
        result = dPolygon(P, vertices)
        np.testing.assert_array_almost_equal(result, expected)

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

        # Expected: positive inside the circle but outside the rectangle,
        # negative inside both shapes, and positive outside the circle
        expected = np.maximum(circle[:, -1], -rectangle[:, -1])

        np.testing.assert_array_almost_equal(result[:, -1], expected)

    def test_dIntersect(self, point_grid):
        circle1 = dCircle(point_grid, -1, 0, 2)
        circle2 = dCircle(point_grid, 1, 0, 2)

        result = dIntersect(circle1, circle2)

        # Expected: negative inside both circles, positive outside either circle
        expected = np.maximum(circle1[:, -1], circle2[:, -1])

        np.testing.assert_array_almost_equal(result[:, -1], expected)

    def test_dUnion(self, point_grid):
        circle = dCircle(point_grid, 2, 2, 1.5)
        rectangle = dRectangle(point_grid, -3, 0, -1, 1)

        result = dUnion(circle, rectangle)

        # Expected: negative inside either shape, positive outside both shapes
        expected = np.minimum(circle[:, -1], rectangle[:, -1])

        np.testing.assert_array_almost_equal(result[:, -1], expected)


if __name__ == "__main__":
    pytest.main()
