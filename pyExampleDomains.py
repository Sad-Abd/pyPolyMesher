import numpy as np
from pyPolyMesher import Domain
from pydFunctions import *


class CookDomain(Domain):
    def __init__(self, name, BdBox, PFix, DistFnc):
        super().__init__(name, BdBox, PFix, DistFnc)

    def compute(self, Demand, BdBox, Arg=None):
        if Demand == 'Dist':
            x = self.DistFnc(Arg, BdBox)
        elif Demand == 'BC':
            x = self.BndryCnds(*Arg, BdBox)
        return x

    def BndryCnds(self, Node, Element, BdBox):
        eps = 0.1 * np.sqrt((BdBox[1] - BdBox[0]) *
                            (BdBox[3] - BdBox[2]) / Node.shape[0])
        LeftsideNodes = np.where(Node[:, 0] < eps)[0]
        Supp = np.ones((LeftsideNodes.shape[0], 3), dtype=int)
        Supp[:, 0] = LeftsideNodes

        RightsideNodes = np.where(Node[:, 0] > 48-eps)[0]
        Load = np.zeros((RightsideNodes.shape[0], 3), dtype=int)
        Load[:, 0] = RightsideNodes
        Load[:, 2] = 20
        return [Supp, Load]


class SuspensionDomain(Domain):
    def __init__(self, name, BdBox, PFix, DistFnc):
        super().__init__(name, BdBox, PFix, DistFnc)

    def compute(self, Demand, BdBox, Arg=None):
        if Demand == 'Dist':
            x = self.DistFnc(Arg, BdBox)
        elif Demand == 'BC':
            x = self.BndryCnds(Arg[0], Arg[1], BdBox)
        return x

    def BndryCnds(self, Node, Element, BdBox):
        CornerCircle = np.sqrt(
            (Node[:, 0] - 2.0) ** 2 + (Node[:, 1] - 2.0) ** 2)
        CornerCircle = np.argsort(CornerCircle)
        UpperCircle = np.sqrt((Node[:, 0] - 2.0) **
                              2 + (Node[:, 1] - 16.0) ** 2)
        UpperCircle = np.argsort(UpperCircle)
        Supp = np.ones((2, 3), dtype=int)
        Supp[0, :] = [CornerCircle[0], 1, 1]
        Supp[1, :] = [UpperCircle[0], 1, 0]
        RightCircle = np.sqrt((Node[:, 0] - 20.0) **
                              2 + (Node[:, 1] - 2.5) ** 2)
        RightCircle = np.argsort(RightCircle)
        Load = np.ones((1, 3))
        Load[0, :] = np.array([RightCircle[0], -8, -1]).reshape((-1, 3))
        x = [Supp, Load]
        return x


class MichellDomain(Domain):
    def __init__(self, name, BdBox, PFix, DistFnc):
        super().__init__(name, BdBox, PFix, DistFnc)

    def compute(self, Demand, BdBox, Arg=None):
        if Demand == 'Dist':
            x = self.DistFnc(Arg, BdBox)
        elif Demand == 'BC':
            x = self.BndryCnds(Arg[0], Arg[1], BdBox)
        return x

    def BndryCnds(self, Node, Element, BdBox):
        eps = 0.1 * ((BdBox[1] - BdBox[0]) *
                     (BdBox[3] - BdBox[2]) / Node.shape[0]) ** 0.5
        CircleNodes = [i for i, (x, y) in enumerate(
            Node) if abs((x ** 2 + y ** 2) ** 0.5 - 1.0) < eps]
        Supp = np.array([[node, 1, 1]
                        for node in CircleNodes]).reshape((-1, 3))
        MidRightFace = [
            ((x - BdBox[1]) ** 2 + (y - (BdBox[2] + BdBox[3]) / 2) ** 2) for x, y in Node]
        MidRightFace = [i for i, _ in sorted(
            enumerate(MidRightFace), key=lambda x: x[1])]
        Load = np.array([MidRightFace[0], 0, -1]).reshape((-1, 3))
        x = [Supp, Load]
        return x


class WrenchDomain(Domain):
    def __init__(self, name, BdBox, PFix, DistFnc):
        super().__init__(name, BdBox, PFix, DistFnc)

    def compute(self, Demand, BdBox, Arg=None):
        if Demand == 'Dist':
            x = self.DistFnc(Arg, BdBox)
        elif Demand == 'BC':
            x = self.BndryCnds(*Arg, BdBox)
        return x

    def BndryCnds(self, Node, Element, BdBox):
        eps = 0.1 * np.sqrt((BdBox[1] - BdBox[0]) *
                            (BdBox[3] - BdBox[2]) / Node.shape[0])
        RightCircleNodes = np.where(
            np.abs(np.sqrt((Node[:, 0] - 2) ** 2 + Node[:, 1] ** 2) - 0.3) < eps)[0]
        Supp = np.ones((RightCircleNodes.shape[0], 3), dtype=int)
        Supp[:, 0] = RightCircleNodes
        LeftHalfCircleNodes = np.where(np.abs(np.maximum(
            np.sqrt(Node[:, 0] ** 2 + Node[:, 1] ** 2) - 0.175, Node[:, 1])) < eps)[0]
        Load = -0.1 * np.ones((LeftHalfCircleNodes.shape[0], 3))
        Load[:, 0] = LeftHalfCircleNodes
        Load[:, 1] = 0
        x = [Supp, np.array(Load)]
        return x


class HornDomain(Domain):
    def __init__(self, name, BdBox, PFix, DistFnc):
        super().__init__(name, BdBox, PFix, DistFnc)

    def compute(self, Demand, BdBox, Arg=None):
        if Demand == 'Dist':
            x = self.DistFnc(Arg, BdBox)
        elif Demand == 'BC':
            x = self.BndryCnds(Arg[0], Arg[1], BdBox)
        return x

    def BndryCnds(self, Node, Element, BdBox):
        x = [None, None]  # No boundary conditions specified for this problem
        return x


class MbbDomain(Domain):
    def __init__(self, name, BdBox, PFix, DistFnc):
        super().__init__(name, BdBox, PFix, DistFnc)

    def compute(self, Demand, BdBox, Arg=None):
        if Demand == "Dist":
            x = self.DistFnc(Arg, BdBox)
        elif Demand == "BC":
            x = self.BndryCnds(*Arg, BdBox)
        return x

    def DistFnc(self, P, BdBox):
        Dist = dRectangle(P, BdBox[0], BdBox[1], BdBox[2], BdBox[3])
        return Dist

    def BndryCnds(self, Node, Element, BdBox):
        eps = 0.1 * np.sqrt((BdBox[1] - BdBox[0]) *
                            (BdBox[3] - BdBox[2]) / Node.shape[0])
        LeftEdgeNodes = np.where(np.abs(Node[:, 0] - BdBox[0]) < eps)[0]
        LeftUpperNode = np.where(
            np.logical_and(
                np.abs(Node[:, 0] - BdBox[0]
                       ) < eps, np.abs(Node[:, 1] - BdBox[3]) < eps
            )
        )[0]
        RightBottomNode = np.where(
            np.logical_and(
                np.abs(Node[:, 0] - BdBox[1]
                       ) < eps, np.abs(Node[:, 1] - BdBox[2]) < eps
            )
        )[0]
        FixedNodes = np.concatenate((LeftEdgeNodes, RightBottomNode))
        Supp = np.zeros((len(FixedNodes), 3), dtype=int)
        Supp[:, 0] = FixedNodes
        Supp[:-1, 1] = 1
        Supp[-1, 2] = 1
        Load = np.zeros((1, 3))
        Load[0, 0], Load[0, 1], Load[0, 2] = LeftUpperNode[0], 0, -0.5
        x = [Supp, Load]
        return x
