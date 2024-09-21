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

    def test_dDiff(self):
        d1 = np.array([[1, 1], [2, 2]])
        d2 = np.array([[1, 1], [2, 2]])
        expected = np.array([[1, 1, 1], [2, 2, 2]])
        result = dDiff(d1, d2)
        np.testing.assert_array_almost_equal(result, expected)

    def test_dIntersect(self):
        d1 = np.array([[1, 1], [2, 2]])
        d2 = np.array([[1, 1], [2, 2]])
        expected = np.array([[1, 1, 1], [2, 2, 2]])
        result = dIntersect(d1, d2)
        np.testing.assert_array_almost_equal(result, expected)

    def test_dUnion(self):
        d1 = np.array([[1, 1], [2, 2]])
        d2 = np.array([[1, 1], [2, 2]])
        expected = np.array([[1, 1, 1], [2, 2, 2]])
        result = dUnion(d1, d2)
        np.testing.assert_array_almost_equal(result, expected)


if __name__ == "__main__":
    pytest.main()
