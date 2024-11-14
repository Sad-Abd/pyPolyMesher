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
    - dLineExact(P, x1, y1, x2, y2): Calculate the exact signed distance 
      from points to a line segment.
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


def dLineExact(P, x1, y1, x2, y2):
    """
    Calculate the exact signed distance from points P to a line segment defined
    by two endpoints (x1, y1) and (x2, y2).

    This function accurately handles cases where points are close to the line segment
    or its endpoints, by projecting points onto the line segment itself and calculating
    distances to the closest endpoint for points outside the segment.

    Parameters:
        P (numpy.ndarray): An array of 2D points (shape: (N, 2)).
        x1 (float): X-coordinate of the first endpoint of the line segment.
        y1 (float): Y-coordinate of the first endpoint of the line segment.
        x2 (float): X-coordinate of the second endpoint of the line segment.
        y2 (float): Y-coordinate of the second endpoint of the line segment.

    Returns:
        numpy.ndarray: An array of signed distances from each point in P to the line segment.
    """
    # Vector from point (x1, y1) to point (x2, y2)
    line_vec = np.array([x2 - x1, y2 - y1])
    line_len = np.linalg.norm(line_vec)
    line_unitvec = line_vec / line_len

    # Vector from point (x1, y1) to each point in P
    vec_to_point = P - np.array([x1, y1])

    # Project the vector to the line and calculate distance
    proj_lengths = np.dot(vec_to_point, line_unitvec)
    proj_points = np.outer(proj_lengths, line_unitvec) + np.array([x1, y1])

    # Calculate the distances to the line segment
    d = np.linalg.norm(proj_points - P, axis=1)

    # Find points that are outside the segment and calculate distances to the closest endpoint
    outside_start = proj_lengths < 0
    outside_end = proj_lengths > line_len
    d[outside_start] = np.linalg.norm(P[outside_start] - np.array([x1, y1]), axis=1)
    d[outside_end] = np.linalg.norm(P[outside_end] - np.array([x2, y2]), axis=1)

    # Calculate signed distances
    perp_vec = np.array([line_unitvec[1], -line_unitvec[0]])
    sign = np.sign(np.dot(vec_to_point, perp_vec))
    signed_distances = d * sign

    return np.column_stack((signed_distances, signed_distances))


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


def dEllipse(P, xc, yc, a, b, theta=0.0):
    """
    Compute the signed distance from points P to an ellipse centered at (xc, yc),
    with semi-axes defined by a and b, and optionally rotated by theta.

    Parameters:
    -----------
    P : numpy.ndarray
        Array of points with shape (N, 2) where each row is a point (x, y).
    xc, yc : float
        Coordinates of the ellipse center.
    a, b : float
        Semi-major and semi-minor axes of the ellipse.
    theta : float, optional
        Rotation angle of the ellipse in radians (default: 0.).

    Returns:
    --------
    numpy.ndarray
        Column stack of the signed distances for each point.
    """

    # Translate and rotate points
    P_translated = P - np.array([xc, yc])

    if theta != 0:
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        rotation_matrix = np.array([[cos_theta, -sin_theta], [sin_theta, cos_theta]])
        P_translated = P_translated @ rotation_matrix.T

    e = np.array([a, b])
    ei = 1.0 / e
    e2 = e * e
    ve = ei * np.array([e2[0] - e2[1], e2[1] - e2[0]])
    distances = np.zeros(P.shape[0])

    # Iterate over each point
    for i, p in enumerate(P_translated):
        p_abs = np.abs(p)

        # Initialize with a point on the ellipse at a 45-degree angle
        t = np.array([np.sqrt(2) / 2, np.sqrt(2) / 2])

        # Iteratively refine the estimate of the closest point on the ellipse
        for _ in range(3):
            v = ve * t * t * t
            u = (p_abs - v) / np.linalg.norm(p_abs - v) * np.linalg.norm(t * e - v)
            w = ei * (v + u)
            t = np.clip(w, 0.0, 1.0)
            t /= np.linalg.norm(t)

        # Calculate the nearest point and signed distance
        nearest_abs = t * e
        dist = np.linalg.norm(p_abs - nearest_abs)

        # Set distance as negative if the point is inside the ellipse
        distances[i] = (
            -dist if np.dot(p_abs, p_abs) < np.dot(nearest_abs, nearest_abs) else dist
        )

    return np.column_stack((distances, distances))


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


def _is_point_in_polygon(P, vertices):
    """
    Determine if points are inside a polygon using the ray-casting algorithm.

    Parameters:
        P (numpy.ndarray): An array of 2D points (shape: (N, 2)).
        vertices (list): A list of vertices defining the polygon.

    Returns:
        numpy.ndarray: An array of boolean values indicating if each point is inside the polygon.
    """
    px, py = P[:, 0], P[:, 1]
    n = len(vertices)

    inside = np.zeros(px.shape, dtype=bool)

    for i in range(n):
        x1, y1 = vertices[i]  # Get the coordinates of the current vertex
        x2, y2 = vertices[
            (i + 1) % n
        ]  # Get the coordinates of the next vertex (wrapping around)

        # Check if the point's y-coordinate is between the y-coordinates of the current edge's vertices
        cond1 = (y1 > py) != (y2 > py)

        # Handle the case where the edge is horizontal
        if y1 == y2:
            # If the point's y-coordinate matches the horizontal edge, check if it's within the x-range of the edge
            cond2 = (py == y1) & (np.minimum(x1, x2) <= px) & (px <= np.maximum(x1, x2))
        else:
            # Calculate the x-coordinate of the intersection of the polygon edge with a horizontal line through the point
            # Then check if the point's x-coordinate is to the left of this intersection
            cond2 = px < (x2 - x1) * (py - y1) / (y2 - y1) + x1

        # If both conditions are met, toggle the 'inside' flag
        inside ^= cond1 & cond2

    return inside


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

    # Calculate distance to each edge
    distances = np.array(
        [
            dLineExact(P, *vertices[i], *vertices[i + 1])[:, -1]
            for i in range(len(vertices) - 1)
        ]
    ).T

    # Get the minimum distance (considering absolute values for accurate distance calculation)
    min_distances = np.min(np.abs(distances), axis=1)

    # Determine if the points are inside or outside the polygon
    inside = _is_point_in_polygon(P, vertices)

    # Apply the correct sign
    signed_distances = min_distances * (inside * -2 + 1)

    return np.column_stack((distances, signed_distances))


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
