import numpy as np
from pydFunctions import *

def MichellDomain(Demand, Arg=None):
    BdBox = [0, 5, -2, 2]
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
    d1 = dRectangle(P, BdBox[0], BdBox[1], BdBox[2], BdBox[3])
    d2 = dCircle(P, 0, 0, BdBox[3] / 2)
    Dist = dDiff(d1, d2)
    return Dist

# ---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
def BndryCnds(Node, Element, BdBox):
    eps = 0.1 * ((BdBox[1] - BdBox[0]) * (BdBox[3] - BdBox[2]) / Node.shape[0]) ** 0.5
    CircleNodes = [i for i, (x, y) in enumerate(Node) if abs((x ** 2 + y ** 2) ** 0.5 - 1.0) < eps]
    Supp = np.array([[node, 1, 1] for node in CircleNodes]).reshape((-1,3))
    MidRightFace = [((x - BdBox[1]) ** 2 + (y - (BdBox[2] + BdBox[3]) / 2) ** 2) for x, y in Node]
    MidRightFace = [i for i, _ in sorted(enumerate(MidRightFace), key=lambda x: x[1])]
    Load = np.array([MidRightFace[0], 0, -1]).reshape((-1,3))
    x = [Supp, Load]
    return x

# ----------------------------------------------------- SPECIFY FIXED POINTS
def FixedPoints(BdBox):
    PFix = [5,0]
    return np.array(PFix).reshape((-1,2))