import numpy as np


def dCircle(P, xc, yc, r):
    d = np.sqrt((P[:, 0] - xc) ** 2 + (P[:, 1] - yc) ** 2) - r
    d = np.column_stack((d, d))
    return d


def dDiff(d1, d2):
    d = np.column_stack((d1[:, :-1], d2[:, :-1]))
    d = np.column_stack((d, np.maximum(d1[:, -1], -d2[:, -1])))
    return d


def dIntersect(d1, d2):
    d = np.column_stack((d1[:, :-1], d2[:, :-1]))
    d = np.column_stack((d, np.maximum(d1[:, -1], d2[:, -1])))
    return d


def dLine(P, x1, y1, x2, y2):
    a = np.array([x2 - x1, y2 - y1])
    a = a / np.linalg.norm(a)
    b = P[:, 0:2] - np.array([x1, y1])
    d = np.dot(b, np.array([-a[1], a[0]]))
    d = np.column_stack((d, d))
    return d


def dRectangle(P, x1, x2, y1, y2):
    d = np.column_stack((x1 - P[:, 0], P[:, 0] - x2, y1 - P[:, 1], P[:, 1] - y2))
    d = np.column_stack((d, np.max(d, axis=1)))
    return d


def dUnion(d1, d2):
    d = np.column_stack((d1[:, :-1], d2[:, :-1]))
    d = np.column_stack((d, np.minimum(d1[:, -1], d2[:, -1])))
    return d
