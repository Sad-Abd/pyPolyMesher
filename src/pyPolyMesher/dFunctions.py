"""
This module provides functions for computing signed distance fields and operations
on them.

A signed distance field is a grid of values where each value represents the signed
distance from a point on the grid to a specified geometric shape or object.
Positive values typically indicate points outside the shape, and negative values
indicate points inside the shape. These functions can be used for various geometric
computations.

Functions:
    - dLine(P, x1, y1, x2, y2): Calculate signed distances from points to a line
      segment.
    - dCircle(P, xc, yc, r): Calculate signed distances from points to a circle.
    - dRectangle(P, x1, x2, y1, y2): Calculate signed distances from points to a
      rectangle.
    - dPolygon(P, vertices): Calculate the signed distance from points P to
      a polygon defined by its vertices.
    - dDiff(d1, d2): Compute the signed distance field resulting from the difference
      of two fields.
    - dIntersect(*ds): Compute the signed distance field resulting from the
      intersection of two or more fields.
    - dUnion(*ds): Compute the signed distance field resulting from the union of
      two or more fields.

These functions are designed for use in various geometric and computational
geometry applications. For more detailed information on each function, refer to their individual
docstrings.



Notes:
- These function names have been chosen to align with corresponding functions in
  MATLAB, both as a mark of respect and loyalty to the authors of the reference code and
  to make the transition from MATLAB to this code as seamless as possible for users.
"""

import numpy as np

# Shapes:


def dLine(P, x1, y1, x2, y2):
    """
    Calculate the signed distance from points P to a line segment defined
    by two endpoints (x1, y1) and (x2, y2).

    Parameters:
        P (numpy.ndarray): An array of 2D points (shape: (N, 2)).
        x1 (float): X-coordinate of the first endpoint of the line segment.
        y1 (float): Y-coordinate of the first endpoint of the line segment.
        x2 (float): X-coordinate of the second endpoint of the line segment.
        y2 (float): Y-coordinate of the second endpoint of the line segment.

    Returns:
        numpy.ndarray: An array of signed distances from each point in P to the line segment.
    """
    a = np.array([x2 - x1, y2 - y1])
    a = a / np.linalg.norm(a)
    b = P - np.array([x1, y1])
    d = np.dot(b, np.array([a[1], -a[0]]))
    d = np.column_stack((d, d))
    return d


def dCircle(P, xc, yc, r):
    """
    Calculate the signed distance from points P to a circle defined by
    its center (xc, yc) and radius (r).

    Parameters:
        P (numpy.ndarray): An array of 2D points (shape: (N, 2)).
        xc (float): X-coordinate of the circle's center.
        yc (float): Y-coordinate of the circle's center.
        r (float): Radius of the circle.

    Returns:
        numpy.ndarray: An array of signed distances from each point in P to the circle.
    """
    d = np.sqrt((P[:, 0] - xc) ** 2 + (P[:, 1] - yc) ** 2) - r
    d = np.column_stack((d, d))
    return d


def dRectangle(P, x1, x2, y1, y2):
    """
    Calculate the signed distance from points P to a rectangle defined
    by its bottom-left (x1, y1) and top-right (x2, y2) coordinates.

    Parameters:
        P (numpy.ndarray): An array of 2D points (shape: (N, 2)).
        x1 (float): X-coordinate of the bottom-left corner of the rectangle.
        x2 (float): X-coordinate of the top-right corner of the rectangle.
        y1 (float): Y-coordinate of the bottom-left corner of the rectangle.
        y2 (float): Y-coordinate of the top-right corner of the rectangle.

    Returns:
        numpy.ndarray: An array of signed distances from each point in P to the rectangle.
    """
    d = np.column_stack((x1 - P[:, 0], P[:, 0] - x2, y1 - P[:, 1], P[:, 1] - y2))
    d = np.column_stack((d, np.max(d, axis=1)))
    return d


def _is_counter_clockwise(points):
    """
    Determines if the given points form a counter-clockwise ordering in a 2D plane.

    Parameters:
        points (list of tuples): A list of (x, y) coordinate tuples representing the polygon's vertices.

    Returns:
        bool: True if the points are arranged in a counter-clockwise order, False otherwise.
    """
    sum_cross_product = 0
    for i in range(len(points)):
        x1, y1 = points[i]
        x2, y2 = points[(i + 1) % len(points)]
        sum_cross_product += (x2 - x1) * (y2 + y1)
    return sum_cross_product > 0


def dPolygon(P, vertices):
    """
    Calculate the signed distance from points P to a polygon defined by its vertices.

    Parameters:
        P (numpy.ndarray): An array of 2D points (shape: (N, 2)).
        vertices (list): A list of vertices defining the polygon.

    Returns:
        numpy.ndarray: An array of signed distances from each point in P to the polygon.
    """
    if _is_counter_clockwise(vertices):
        vertices = vertices[::-1]
    if vertices[0] != vertices[-1]:
        vertices.append(vertices[0])
    d = np.column_stack(
        [
            dLine(P, *vertices[i], *vertices[i + 1])[:, -1]
            for i in range(len(vertices) - 1)
        ]
    )
    d = np.column_stack((d, np.max(d, axis=1)))
    return d


# Boolean Operations:


def dDiff(d1, d2):
    """
    Calculate the signed distance field resulting from the difference of
    two distance fields (d1 and d2).

    Parameters:
        d1 (numpy.ndarray): The first distance field.
        d2 (numpy.ndarray): The second distance field.

    Returns:
        numpy.ndarray: The resulting signed distance field.
    """
    d = np.column_stack((d1[:, :-1], d2[:, :-1]))
    d = np.column_stack((d, np.maximum(d1[:, -1], -d2[:, -1])))
    return d


def dIntersect(*ds):
    """
    Calculate the signed distance field resulting from the intersection
    of two or more distance fields.

    Parameters:
        *ds (numpy.ndarray): Distance fields to be intersected.

    Returns:
        numpy.ndarray: The resulting signed distance field.
    """

    d = np.column_stack([each[:, :-1] for each in ds])
    d = np.column_stack((d, np.maximum.reduce([each[:, -1] for each in ds])))

    return d


def dUnion(*ds):
    """
    Calculate the signed distance field resulting from the union
    of two or more distance fields.

    Parameters:
        *ds (numpy.ndarray): Distance fields to be intersected.

    Returns:
        numpy.ndarray: The resulting signed distance field.
    """
    d = np.column_stack([each[:, :-1] for each in ds])
    d = np.column_stack((d, np.minimum.reduce([each[:, -1] for each in ds])))
    return d
