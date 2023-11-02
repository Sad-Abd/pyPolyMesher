import numpy as np
from pyPolyMesher import *
from pyExampleDomains import *
from pydFunctions import *


# for structured mesh
# nelx = 5
# nely = 4
# dx = 3 / nelx
# dy = 1 / nely
# x_range = np.arange(dx / 2, 3, dx)
# y_range = np.arange(dy / 2, 1, dy)
# X, Y = np.meshgrid(x_range, y_range)
# P = np.column_stack((X.ravel(), Y.ravel()))
# PolyMesher(MbbDomain, 20, 30, P)

# for unstructured mesh with given points (can be used to compare results with MATLAB code)
# P = np.array([[2.60026522, 0.03321281],
#  [0.01350291, 0.39790721],
#  [1.06513347, 0.40730687],
#  [2.54716707, 0.76359476],
#  [1.29260857, 0.39074023],
#  [2.41308243, 0.38506563],
#  [2.22468734, 0.69673122],
#  [0.95129956, 0.09203691],
#  [0.6433727,  0.67682601],
#  [2.1696533,  0.3882989 ]])
# PolyMesher(MbbDomain, 10, 50, P)

# for unstructured mesh

# ---------------------------------------------------------------------------------

# def DistFnc(P, BdBox):
#     Dist = dRectangle(P, BdBox[0], BdBox[1], BdBox[2], BdBox[3])
#     return Dist


# BdBox = [0, 3, 0, 1]
# PFix = []

# mbb_domain = MbbDomain("Mbb Domain", BdBox, PFix, DistFnc)
# PolyMesher(mbb_domain, 50, 100)

# ---------------------------------------------------------------------------------

# def DistFnc(P, BdBox):
#     d1 = dCircle(P, 0, 0, 1)
#     d2 = dCircle(P, -0.4, 0, 0.55)
#     d3 = dLine(P, 0, 0, 1, 0)
#     Dist = dIntersect(d3, dDiff(d1, d2))
#     return Dist


# BdBox = [-1, 1, 0, 1]
# PFix = []

# horn_domain = HornDomain("Horn Domain", BdBox, PFix, DistFnc)
# PolyMesher(horn_domain, 150, 50, anim=True)

# ---------------------------------------------------------------------------------

# def DistFnc(P, BdBox):
#     d1 = dLine(P, 0, 0.3, 0, -0.3)
#     d2 = dLine(P, 0, -0.3, 2, -0.5)
#     d3 = dLine(P, 2, -0.5, 2, 0.5)
#     d4 = dLine(P, 2, 0.5, 0, 0.3)
#     d5 = dCircle(P, 0, 0, 0.3)
#     d6 = dCircle(P, 2, 0, 0.5)
#     douter = dUnion(d6, dUnion(d5,
#                                dIntersect(d4, dIntersect(d3, dIntersect(d2, d1)))))
#     d7 = dCircle(P, 0, 0, 0.175)
#     d8 = dCircle(P, 2, 0, 0.3)
#     din = dUnion(d8, d7)
#     Dist = dDiff(douter, din)
#     return Dist


# BdBox = [-0.3, 2.5, -0.5, 0.5]
# PFix = []

# wrench_domain = WrenchDomain("Wrench Domain", BdBox, PFix, DistFnc)
# PolyMesher(wrench_domain, 150, 100)

# ---------------------------------------------------------------------------------

# def DistFnc(P, BdBox):
#     d1 = dRectangle(P, BdBox[0], BdBox[1], BdBox[2], BdBox[3])
#     d2 = dCircle(P, 0, 0, BdBox[3] / 2)
#     Dist = dDiff(d1, d2)
#     return Dist


# BdBox = [0, 5, -2, 2]
# PFix = [5, 0]

# michell_domain = MichellDomain("Michell Domain", BdBox, PFix, DistFnc)
# PolyMesher(michell_domain, 20, 100)

# ---------------------------------------------------------------------------------

# def DistFnc(P, BdBox):
#     d1 = dRectangle(P, 0, 18.885, 0, 14.56)
#     d2 = dLine(P, 18.885, 1.3030, 4, 0)
#     d3 = dLine(P, 3.92, 14.56, 6.1699, 6.88)
#     d4 = dLine(P, 9.8651, 4.0023, 18.885, 3.70)
#     d5 = dLine(P, 4, 0, 0, 4)
#     d13 = dLine(P, 0, 14, 3.92, 14.56)
#     d14 = dCircle(P, 10, 8, 4)
#     d15 = dLine(P, 9.8651, 4.0023, 6.1699, 6.88)
#     d = dDiff(dDiff(dDiff(dDiff(d1, d2), d5), d13),
#               dUnion(dDiff(dIntersect(d3, d4), d15), d14))
#     d6 = dCircle(P, 2, 2, 2)
#     d7 = dCircle(P, 4, 2, 2)
#     d8 = dCircle(P, 2, 4, 2)
#     d = dUnion(d, dUnion(d6, dUnion(d7, d8)))
#     d9 = dCircle(P, 2, 14, 2)
#     d10 = dCircle(P, 2, 16, 2)
#     d11 = dCircle(P, 18.885, 2.5, 1.2)
#     d12 = dCircle(P, 20, 2.5, 1.2)
#     Dist = dUnion(d, dUnion(d9, dUnion(d10, dUnion(d11, d12))))
#     return Dist


# BdBox = [-2, 24, -2, 24]
# PFix = [[2, 2], [2, 16], [20, 2.5]]

# suspension_domain = SuspensionDomain("Suspension Domain", BdBox, PFix, DistFnc)
# PolyMesher(suspension_domain, 750, 150)

# ---------------------------------------------------------------------------------

def DistFnc(P, BdBox):
    d1 = dLine(P, 0., 44., 0., 0.)
    d2 = dLine(P, 0., 0., 48., 44.)
    d3 = dLine(P, 48., 44., 48., 60.)
    d4 = dLine(P, 48., 60., 0., 44.)
    Dist = dIntersect(d4, dIntersect(d3, dIntersect(d2, d1)))
    return Dist


BdBox = [0, 48, 0, 60]
PFix = []

cook_domain = CookDomain("Cook Domain", BdBox, PFix, DistFnc)
Node, Element, Supp, Load, P = PolyMesher(cook_domain, 50, 20)
