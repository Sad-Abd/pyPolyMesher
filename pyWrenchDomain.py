import numpy as np
from pydFunctions import *

def WrenchDomain(Demand, Arg=None):
    BdBox = [-0.3, 2.5, -0.5, 0.5]
    if Demand == 'Dist':
        x = DistFnc(Arg, BdBox)
    elif Demand == 'BC':
        x = BndryCnds(*Arg, BdBox)
    elif Demand == 'BdBox':
        x = BdBox
    elif Demand == 'PFix':
        x = FixedPoints(BdBox)
    return x

def DistFnc(P, BdBox):
    d1 = dLine(P, 0, 0.3, 0, -0.3)
    d2 = dLine(P, 0, -0.3, 2, -0.5)
    d3 = dLine(P, 2, -0.5, 2, 0.5)
    d4 = dLine(P, 2, 0.5, 0, 0.3)
    d5 = dCircle(P, 0, 0, 0.3)
    d6 = dCircle(P, 2, 0, 0.5)
    douter = dUnion(d6, dUnion(d5,
               dIntersect(d4, dIntersect(d3, dIntersect(d2, d1)))))
    d7 = dCircle(P, 0, 0, 0.175)
    d8 = dCircle(P, 2, 0, 0.3)
    din = dUnion(d8, d7)
    Dist = dDiff(douter, din)
    return Dist

def BndryCnds(Node, Element, BdBox):
    eps = 0.1 * np.sqrt((BdBox[1] - BdBox[0]) * (BdBox[3] - BdBox[2]) / Node.shape[0])
    RightCircleNodes = np.where(np.abs(np.sqrt((Node[:, 0] - 2) ** 2 + Node[:, 1] ** 2) - 0.3) < eps)[0]
    Supp = np.ones((RightCircleNodes.shape[0], 3), dtype=int)
    Supp[:, 0] = RightCircleNodes
    LeftHalfCircleNodes = np.where(np.abs(np.maximum(np.sqrt(Node[:, 0] ** 2 + Node[:, 1] ** 2) - 0.175, Node[:, 1])) < eps)[0]
    Load = -0.1 * np.ones((LeftHalfCircleNodes.shape[0], 3))
    Load[:, 0] = LeftHalfCircleNodes
    Load[:, 1] = 0
    x = [Supp, np.array(Load)]
    return x

def FixedPoints(BdBox):
    PFix = []
    return np.array(PFix)
