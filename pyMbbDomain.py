import numpy as np
from pydFunctions import dRectangle


def MbbDomain(Demand, Arg=None):
    BdBox = [0, 3, 0, 1]
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
    Dist = dRectangle(P, BdBox[0], BdBox[1], BdBox[2], BdBox[3])
    return Dist


def BndryCnds(Node, Element, BdBox):
    eps = 0.1 * np.sqrt((BdBox[1] - BdBox[0]) * (BdBox[3] - BdBox[2]) / Node.shape[0])
    LeftEdgeNodes = np.where(np.abs(Node[:, 0] - BdBox[0]) < eps)[0]
    LeftUpperNode = np.where(
        np.logical_and(
            np.abs(Node[:, 0] - BdBox[0]) < eps, np.abs(Node[:, 1] - BdBox[3]) < eps
        )
    )[0]
    RightBottomNode = np.where(
        np.logical_and(
            np.abs(Node[:, 0] - BdBox[1]) < eps, np.abs(Node[:, 1] - BdBox[2]) < eps
        )
    )[0]
    FixedNodes = np.concatenate((LeftEdgeNodes, RightBottomNode))
    Supp = np.zeros((len(FixedNodes), 3), dtype=int)
    Supp[:, 0] = FixedNodes
    Supp[:-1, 1] = 1
    Supp[-1, 2] = 1
    Load = np.zeros((1, 3))
    Load[0,0], Load[0,1], Load[0,2] = LeftUpperNode[0], 0, -0.5 
    x = [Supp, Load]
    return x


def FixedPoints(BdBox):
    PFix = []
    return np.array(PFix)
