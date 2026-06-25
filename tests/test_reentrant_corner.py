import numpy as np

from pyPolyMesher import reentrant_corner_seeds


def _reflect_line(p, corner, direction):
    """Reflect point p about the line through `corner` with unit `direction`."""
    corner = np.asarray(corner, float)
    direction = np.asarray(direction, float)
    direction = direction / np.linalg.norm(direction)
    v = np.asarray(p, float) - corner
    return corner + 2.0 * np.dot(v, direction) * direction - v


def test_returns_three_seeds():
    seeds = reentrant_corner_seeds([1, 1], [1, 0], [0, 1], 0.1)
    assert len(seeds) == 3
    assert all(len(s) == 2 for s in seeds)


def test_seeds_equidistant_from_corner():
    corner = np.array([1.0, 1.0])
    seeds = np.array(reentrant_corner_seeds(corner, [1, 0], [0, 1], 0.12))
    dists = np.linalg.norm(seeds - corner, axis=1)
    assert np.allclose(dists, dists[0])


def test_all_seeds_share_common_external_reflection():
    # Paper (Fig. 15): y1, y2, y3 are fixed such that they all reflect onto the
    # same external point y in the notch.
    corner = np.array([1.0, 1.0])
    e1, e2, delta = [1, 0], [0, 1], 0.15
    y = corner + delta * (np.array(e1) + np.array(e2))
    y1, y2, y3 = (np.array(s) for s in reentrant_corner_seeds(corner, e1, e2, delta))

    assert np.allclose(_reflect_line(y1, corner, e1), y)   # reflect across wall 1
    assert np.allclose(_reflect_line(y2, corner, e2), y)   # reflect across wall 2
    assert np.allclose(2.0 * corner - y3, y)               # reflect through corner


def test_seeds_lie_inside_lshape_domain():
    # For the L-shape (notch = top-right square), all three seeds must be inside.
    seeds = np.array(reentrant_corner_seeds([1, 1], [1, 0], [0, 1], 0.12))
    in_notch = (seeds[:, 0] > 1) & (seeds[:, 1] > 1)
    assert not in_notch.any()


def test_direction_inputs_are_normalized():
    # Non-unit direction vectors should give the same result as unit ones.
    a = np.array(reentrant_corner_seeds([0, 0], [3, 0], [0, 5], 0.1))
    b = np.array(reentrant_corner_seeds([0, 0], [1, 0], [0, 1], 0.1))
    assert np.allclose(a, b)
