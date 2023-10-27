# Cook's membrane geometry definition

import numpy as np
from pydFunctions import dIntersect,dLine


def CookDomain(Demand, Arg=None):
    BdBox = [0, 48, 0, 60]
    if Demand == "Dist":
        x = DistFnc(Arg, BdBox)
    elif Demand == "BC":
        x = BndryCnds(*Arg, BdBox)
    elif Demand == "BdBox":
        x = BdBox
    elif Demand == "PFix":
        x = FixedPoints(BdBox)
    return x


def DistFnc(P, BdBox):
    d1 = dLine(P, 0.,44.,0.,0.)
    d2 = dLine(P, 0., 0., 48., 44.)
    d3 = dLine(P, 48., 44., 48., 60.)
    d4 = dLine(P, 48., 60., 0., 44.)
    Dist = dIntersect(d4, dIntersect(d3, dIntersect(d2, d1)))
    return Dist


def BndryCnds(Node, Element, BdBox):
    eps = 0.1 * np.sqrt((BdBox[1] - BdBox[0]) * (BdBox[3] - BdBox[2]) / Node.shape[0])
    LeftsideNodes = np.where(Node[:,0] < eps)[0]
    Supp = np.ones((LeftsideNodes.shape[0], 3), dtype=int)
    Supp[:, 0] = LeftsideNodes

    RightsideNodes = np.where(Node[:,0] > 48-eps)[0]
    Load = np.zeros((RightsideNodes.shape[0], 3), dtype=int)
    Load[:, 0] = RightsideNodes
    Load[:, 2] = 20
    return [Supp, Load]


def FixedPoints(BdBox):
    PFix = []
    return np.array(PFix)