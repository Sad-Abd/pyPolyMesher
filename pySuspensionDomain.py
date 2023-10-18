import numpy as np
from pydFunctions import *

def SuspensionDomain(Demand, Arg=None):
    BdBox = [-2, 24, -2, 24]
    if Demand == 'Dist':
        x = DistFnc(Arg, BdBox)
    elif Demand == 'BC':
        x = BndryCnds(Arg[0], Arg[1], BdBox)
    elif Demand == 'BdBox':
        x = BdBox
    elif Demand == 'PFix':
        x = FixedPoints(BdBox)
    return x

# ----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
def DistFnc(P, BdBox):
    d1 = dRectangle(P, 0, 18.885, 0, 14.56)
    d2 = dLine(P, 18.885, 1.3030, 4, 0)
    d3 = dLine(P, 3.92, 14.56, 6.1699, 6.88)
    d4 = dLine(P, 9.8651, 4.0023, 18.885, 3.70)
    d5 = dLine(P, 4, 0, 0, 4)
    d13 = dLine(P, 0, 14, 3.92, 14.56)
    d14 = dCircle(P, 10, 8, 4)
    d15 = dLine(P, 9.8651, 4.0023, 6.1699, 6.88)
    d = dDiff(dDiff(dDiff(dDiff(d1, d2), d5), d13), dUnion(dDiff(dIntersect(d3, d4), d15), d14))
    d6 = dCircle(P, 2, 2, 2)
    d7 = dCircle(P, 4, 2, 2)
    d8 = dCircle(P, 2, 4, 2)
    d = dUnion(d, dUnion(d6, dUnion(d7, d8)))
    d9 = dCircle(P, 2, 14, 2)
    d10 = dCircle(P, 2, 16, 2)
    d11 = dCircle(P, 18.885, 2.5, 1.2)
    d12 = dCircle(P, 20, 2.5, 1.2)
    Dist = dUnion(d, dUnion(d9, dUnion(d10, dUnion(d11, d12))))
    return Dist

# ---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
def BndryCnds(Node, Element, BdBox):
    CornerCircle = np.sqrt((Node[:, 0] - 2.0) ** 2 + (Node[:, 1] - 2.0) ** 2)
    CornerCircle = np.argsort(CornerCircle)
    UpperCircle = np.sqrt((Node[:, 0] - 2.0) ** 2 + (Node[:, 1] - 16.0) ** 2)
    UpperCircle = np.argsort(UpperCircle)
    Supp = np.ones((2, 3), dtype=int)
    Supp[0, :] = [CornerCircle[0], 1, 1]
    Supp[1, :] = [UpperCircle[0], 1, 0]
    RightCircle = np.sqrt((Node[:, 0] - 20.0) ** 2 + (Node[:, 1] - 2.5) ** 2)
    RightCircle = np.argsort(RightCircle)
    Load = np.ones((1, 3))
    Load[0, :] = np.array([RightCircle[0], -8, -1]).reshape((-1,3))
    x = [Supp, Load]
    return x

# ----------------------------------------------------- SPECIFY FIXED POINTS
def FixedPoints(BdBox):
    PFix = np.array([[2, 2], [2, 16], [20, 2.5]])
    return PFix

