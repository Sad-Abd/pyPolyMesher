import numpy as np

import matplotlib.pyplot as plt

from scipy.spatial import Voronoi, voronoi_plot_2d
import scipy.sparse as sp


def PolyMesher(Domain, NElem, MaxIter, P=None):
    if P is None:
        P = PolyMshr_RndPtSet(NElem, Domain)

    NElem = P.shape[0]
    Tol = 5e-6
    It = 0
    Err = 1
    c = 1.5
    BdBox = Domain("BdBox")
    PFix = Domain("PFix")
    Area = (BdBox[1] - BdBox[0]) * (BdBox[3] - BdBox[2])
    Pc = P.copy()

    while It <= MaxIter and Err > Tol:
        Alpha = c * np.sqrt(Area / NElem)
        P = Pc.copy()
        R_P = PolyMshr_Rflct(P, NElem, Domain, Alpha)
        P, R_P = PolyMshr_FixedPoints(P, R_P, PFix)

        vor = Voronoi(np.vstack((P, R_P)))
        # fig = voronoi_plot_2d(vor)
        # plt.show()
        Node = vor.vertices

        # Reordering regions based on points
        Element = [vor.regions[reg] for reg in vor.point_region]

        Pc, A = PolyMshr_CntrdPly(Element, Node, NElem)
        Area = np.sum(np.abs(A))
        Err = (
            np.sqrt(np.sum((A**2) * np.sum((Pc - P) * (Pc - P), axis=1)))
            * NElem
            / (Area**1.5)
        )
        print(f"It: {It}   Error: {Err}")
        It += 1

        if NElem <= 2000:
            PolyMshr_PlotMsh(Node, Element, NElem)

    Node, Element = PolyMshr_ExtrNds(Node, Element[:NElem])

    Node, Element = PolyMshr_CllpsEdgs(Node, Element, 0.1)

    # to be fixed:
    # Node, Element = PolyMshr_RsqsNds(Node, Element, NElem)

    BC = Domain('BC',[Node, Element])
    Supp = BC[0]
    Load = BC[1]

    PolyMshr_PlotMsh(Node, Element, NElem, wait=True)
    
    return Node, Element, Supp, Load, P


def PolyMshr_RndPtSet(NElem, Domain):
    P = np.zeros((NElem, 2))
    BdBox = Domain("BdBox")
    Ctr = 0

    while Ctr < NElem:
        Y = np.random.rand(NElem, 2)
        Y[:, 0] = (BdBox[1] - BdBox[0]) * Y[:, 0] + BdBox[0]
        Y[:, 1] = (BdBox[3] - BdBox[2]) * Y[:, 1] + BdBox[2]
        d = Domain("Dist", Y)
        I = np.where(d[:, -1] < 0)[0]
        NumAdded = min(NElem - Ctr, len(I))
        P[Ctr : Ctr + NumAdded, :] = Y[I[0:NumAdded], :]
        Ctr += NumAdded

    return P


def PolyMshr_FixedPoints(P, R_P, PFix):
    PP = np.vstack((P, R_P))
    for i in range(PFix.shape[0]):
        B = np.argsort(
            np.sqrt((PP[:, 0] - PFix[i, 0]) ** 2 + (PP[:, 1] - PFix[i, 1]) ** 2)
        )
        for j in range(1, 4):
            n = PP[B[j], :] - PFix[i, :]
            n = n / np.linalg.norm(n)
            PP[B[j], :] = PP[B[j], :] - n * (B[j] - B[0])

    P = PP[0 : P.shape[0], :]
    R_P = PP[P.shape[0] :, :]
    return P, R_P


def PolyMshr_Rflct(P, NElem, Domain, Alpha):
    eps = 1e-8
    eta = 0.9
    d = Domain("Dist", P)
    NBdrySegs = d.shape[1] - 1
    n1 = ((Domain("Dist", P + np.array([[eps, 0]] * NElem))) - d) / eps
    n2 = ((Domain("Dist", P + np.array([[0, eps]] * NElem))) - d) / eps
    I = np.abs(d[:, 0:NBdrySegs]) < Alpha
    P1 = np.broadcast_to(P[:, 0][:, np.newaxis], (P[:, 0].shape[0], NBdrySegs))
    P2 = np.broadcast_to(P[:, 1][:, np.newaxis], (P[:, 1].shape[0], NBdrySegs))

    P1 = P1[I] - 2 * n1[:, 0:NBdrySegs][I] * d[:, 0:NBdrySegs][I]
    P2 = P2[I] - 2 * n2[:, 0:NBdrySegs][I] * d[:, 0:NBdrySegs][I]
    R_P = np.vstack((P1, P2)).T
    d_R_P = Domain("Dist", R_P)
    J = (np.abs(d_R_P[:, -1]) >= eta * np.abs(d[:, 0:NBdrySegs][I])) & (
        d_R_P[:, -1] > 0
    )

    R_P = np.unique(R_P[J, :], axis=0)
    return R_P


def PolyMshr_CntrdPly(Element, Node, NElem):
    centroids = []
    areas = []
    counter = 0
    for vertices in Element:
        if counter >= NElem:
            break
        if -1 in vertices:
            continue
        if len(vertices) >= 3:
            # Extract the vertex coordinates
            polygon_vertices = Node[vertices]

            # Compute the centroid of the polygon
            centroid = np.mean(polygon_vertices, axis=0)
            centroids.append(centroid)

            # Compute the area of the polygon using the shoelace formula
            x = polygon_vertices[:, 0]
            y = polygon_vertices[:, 1]
            area = 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))
            areas.append(area)
            counter += 1

    return np.array(centroids), np.array(areas)


def PolyMshr_ExtrNds(Node0, Element0):
    unique_nodes = np.unique(np.concatenate(Element0))
    cNode = np.arange(len(Node0))
    cNode[~np.in1d(cNode, unique_nodes)] = np.max(unique_nodes)
    Node, Element = PolyMshr_RbldLists(Node0, Element0, cNode)
    return Node, Element


def PolyMshr_CllpsEdgs(Node0, Element0, Tol):
    while True:
        cEdge = []
        for ele in Element0:
            if len(ele) < 4:
                continue  # Cannot collapse triangles
            vx = Node0[ele, 0]
            vy = Node0[ele, 1]
            nv = len(vx)
            beta = np.arctan2(vy - np.sum(vy) / nv, vx - np.sum(vx) / nv)
            beta = np.mod(
                beta[np.roll(np.arange(len(beta)), shift=-1)] - beta, 2 * np.pi
            )
            beta_ideal = 2 * np.pi / len(ele)
            Edge = np.column_stack((ele, np.roll(ele, shift=-1)))
            cEdge.extend(Edge[beta < Tol * beta_ideal, :])

        if len(cEdge) == 0:
            break

        cEdge = np.unique(np.sort(cEdge, axis=1), axis=0)
        cNode = np.arange(len(Node0))
        for i in range(cEdge.shape[0]):
            cNode[cEdge[i, 1]] = cNode[cEdge[i, 0]]
        Node0, Element0 = PolyMshr_RbldLists(Node0, Element0, cNode)
    return Node0, Element0


def PolyMshr_RsqsNds(Node0, Element0, NElem):
    NNode0 = Node0.shape[0]
    NElem0 = len(Element0)

    ElemLnght = [len(e) for e in Element0]
    nn = np.sum(np.array(ElemLnght) ** 2)

    i = np.zeros(nn, dtype=int)
    j = np.zeros(nn, dtype=int)
    s = np.ones(nn)
    index = 0

    for el in range(NElem0):
        eNode = Element0[el]
        ElemSet = np.arange(index, index + ElemLnght[el] ** 2)
        i[ElemSet] = np.kron(eNode, np.ones(ElemLnght[el], dtype=int))
        j[ElemSet] = np.kron(eNode, np.ones(ElemLnght[el], dtype=int))
        index += ElemLnght[el] ** 2

    K = sp.csr_matrix((s, (i, j)), shape=(NNode0, NNode0))
    p = sp.csgraph.reverse_cuthill_mckee(K)

    cNode = np.arange(0, NNode0)
    cNode[p[:NNode0]] = np.arange(0, NNode0)

    Node, Element = PolyMshr_RbldLists(Node0, Element0, cNode)

    return Node, Element


def PolyMshr_RbldLists(Node0, Element0, cNode):
    Element = [None] * len(Element0)
    _, ix, jx = np.unique(cNode, return_index=True, return_inverse=True)

    if not np.array_equal(jx.shape, cNode.shape):
        jx = jx.T

    if Node0.shape[0] > len(ix):
        ix[-1] = max(cNode)

    Node = Node0[ix]

    for ind, ele in enumerate(Element0):
        Element[ind] = np.unique(jx[ele])
        vx = Node[Element[ind], 0]
        vy = Node[Element[ind], 1]
        nv = len(vx)
        angles = np.arctan2(vy - np.sum(vy) / nv, vx - np.sum(vx) / nv)
        iix = np.argsort(angles)
        Element[ind] = Element[ind][iix]

    return Node, Element


def PolyMshr_PlotMsh(Node, Element, NElem, Supp=None, Load=None, wait=False):
    plt.clf()

    Element = Element[:NElem]
    Node_set = set()
    for polygon in Element:
        vx = [Node[i, 0] for i in polygon]
        vy = [Node[i, 1] for i in polygon]
        Node_set.update(polygon)
        plt.fill(vx, vy, "-w", edgecolor="black")

    Node_set = Node[list(Node_set)]
    plt.plot(Node_set[:, 0], Node_set[:, 1], "bo", markersize=8)

    if Supp is not None and len(Supp) > 0:  # Plot Supp BC if specified
        plt.plot(Node[Supp[:, 0], 0], Node[Supp[:, 0], 1], "b>", markersize=8)

    if Load is not None and len(Load) > 0:  # Plot Load BC if specified
        plt.plot(Node[Load[:, 0], 0], Node[Load[:, 0], 1], "m^", markersize=8)

    if not wait:
        plt.show(block=False)
        plt.pause(3e-1)
    else:
        plt.show()
