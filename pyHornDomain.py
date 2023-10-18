import numpy as np
from pydFunctions import *

def HornDomain(Demand, Arg=None):
    BdBox = [-1, 1, 0, 1]
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
    d1 = dCircle(P, 0, 0, 1)
    d2 = dCircle(P, -0.4, 0, 0.55)
    d3 = dLine(P, 0, 0, 1, 0)
    Dist = dIntersect(d3, dDiff(d1, d2))
    return Dist

# ---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
def BndryCnds(Node, Element, BdBox):
    x = [None, None]  # No boundary conditions specified for this problem
    return x

# ----------------------------------------------------- SPECIFY FIXED POINTS
def FixedPoints(BdBox):
    PFix = []
    return np.array(PFix)