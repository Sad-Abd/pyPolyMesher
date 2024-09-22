"""
PolyMesher Module

This module provides functions for generating polygon meshes using Lloyd's algorithm. 
The generated polygon meshes are suitable for various computational geometry and finite element analysis applications.

Functions:
    - PolyMesher(Domain, NElem, MaxIter, P=None, anim=False): Perform polygon mesh generation using Lloyd's algorithm.
    - PolyMshr_RndPtSet(NElem, Domain): Generate an initial random point set of size 'NElem' for polygon mesh generation.
    - PolyMshr_FixedPoints(P, R_P, PFix): Adjust points based on fixed points to maintain mesh quality.
    - PolyMshr_Rflct(P, NElem, Domain, Alpha): Reflect points at the boundary for mesh generation.
    - PolyMshr_CntrdPly(Element, Node, NElem): Compute centroids and areas for elements in the mesh.
    - PolyMshr_ExtrNds(Node0, Element0): Extract unique nodes and rebuild node and element lists.
    - PolyMshr_CllpsEdgs(Node0, Element0, Tol): Collapse small edges based on a specified tolerance.
    - PolyMshr_RsqsNds(Node0, Element0, NElem): Rearrange nodes to improve mesh quality using (RCM) algorithm.
    - PolyMshr_RbldLists(Node0, Element0, cNode): Rebuild node and element lists based on node mapping.
    - PolyMshr_PlotMsh(Node, Element, NElem, Supp=None, Load=None, wait=False): Plot the polygon mesh with optional boundary conditions.

Classes:
    - Domain: Represents a mathematically defined domain for polygonal mesh. See the class docstring for details.

These functions and the 'Domain' class can be used for polygon mesh generation in computational 
geometry applications. For more detailed information on each function and the 'Domain' class, 
refer to their individual docstrings.

Notes:
- These function names have been chosen to align with corresponding functions in
  MATLAB, both as a mark of respect and loyalty to the authors of the reference code and
  to make the transition from MATLAB to this code as seamless as possible for users.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.sparse import csr_matrix, csgraph

from . import progress


class Domain:
    """
    Represents a mathematically defined domain for polygonal mesh.

    This class defines a mathematically defined domain based on its name, bounding box,
    signed distance function (SDF), boundary conditions (BC), and fixed points (PFix).

    Attributes:
        name (str): The name of the mathematically defined domain.
        BdBox (list): The bounding box of the domain, defined as [xmin, xmax, ymin, ymax].
        SDF (callable): The signed distance function that provides the distance values.
        BC (callable, optional): The function for setting boundary conditions. Default is None.
        PFix (list, optional): List of fixed points within the domain. Default is an empty list.

    Methods:
        DistFnc(P):
            Computes the distance value for a point P using the signed distance function.

        BndryCnds(Node):
            Determines the boundary conditions for a given node within the domain.

        Plot(n=1000):
            Plots the domain based on the signed distance function.

        CalculateArea(n=1_000_000):
            Calculates the approximate area of the domain using the Monte Carlo method
    """

    def __init__(self, name, BdBox, SDF, BC=None, PFix=[]):
        """
        Initializes a Domain object with the provided attributes.

        Args:
            name (str): The name of the computational domain.
            BdBox (list): The bounding box of the domain, defined as [xmin, xmax, ymin, ymax].
            SDF (callable): The signed distance function that provides the distance values.
            BC (callable, optional): The function for setting boundary conditions. Default is None.
            PFix (list, optional): List of fixed points within the domain. Default is an empty list.
        """
        self.name = name
        self.BdBox = BdBox
        self.PFix = PFix
        self.SDF = SDF
        self.BC = BC

    def DistFnc(self, P):
        """
        Computes the distance value for a point P using the signed distance function.

        Args:
            P (numpy.ndarray): A point or an array of points for which the distance is to be calculated.

        Returns:
            numpy.ndarray: An array of distance values corresponding to the input points.
        """
        return self.SDF(P)

    def BndryCnds(self, Node):
        """
        Determines the boundary conditions for a given node within the domain.

        Args:
            Node (tuple): The coordinates of a node (x, y).

        Returns:
            list: A list containing boundary conditions for the node. It may include None values.
        """
        if self.BC == None:
            return [None, None]
        return self.BC(Node, self.BdBox)

    def Plot(self, n=1000):
        """
        Plots the domain based on the signed distance function.

        Args:
            n (int, optional): The number of points for plotting. Default is 1000.

        Displays:
            A plot of the domain based on the signed distance function.
        """
        x, y = np.meshgrid(
            np.linspace(self.BdBox[0], self.BdBox[1], n),
            np.linspace(self.BdBox[2], self.BdBox[3], n),
        )
        points = np.hstack([x.reshape((-1, 1)), y.reshape((-1, 1))])
        sdf = self.SDF(points)[:, -1]

        inner = np.where(sdf <= 0, 1, 0)

        _, ax = plt.subplots(figsize=(8, 6))
        _ = ax.imshow(
            inner.reshape((n, n)),
            extent=(self.BdBox[0], self.BdBox[1], self.BdBox[2], self.BdBox[3]),
            origin="lower",
            cmap="Purples",
            alpha=0.8,
        )
        ax.contour(x, y, sdf.reshape((n, n)), levels=[0], colors="gold", linewidths=2)
        ax.set_xlabel("X", fontweight="bold")
        ax.set_ylabel("Y", fontweight="bold")
        ax.set_title("Domain Visualization", fontweight="bold", fontsize=16)
        ax.set_aspect("equal")

        plt.show()

    def CalculateArea(self, n=1_000_000):
        """
        Calculates the approximate area of the domain using the Monte Carlo method.

        Args:
            n (int, optional): The number of random points to use. Default is 1,000,000.

        Returns:
            float: The calculated approximate area of the domain.
        """
        xmin, xmax, ymin, ymax = self.BdBox
        total_area = (xmax - xmin) * (ymax - ymin)

        # Generate random points within the bounding box
        points = np.random.uniform(low=[xmin, ymin], high=[xmax, ymax], size=(n, 2))

        # Evaluate SDF for all points at once
        sdf_values = self.SDF(points)[:, -1]

        # Count points inside the domain
        points_inside = np.sum(sdf_values <= 0)

        # Calculate the ratio and multiply by total area
        area_ratio = points_inside / n
        approximate_area = area_ratio * total_area

        return approximate_area


def PolyMesher(Domain, NElem, MaxIter, P=None, anim=False):
    """
    Perform polygon mesh generation using Lloyd's algorithm.

    Args:
        Domain (function): The domain in which the mesh is generated.
        NElem (int): The desired number of elements in the mesh.
        MaxIter (int): The maximum number of iterations for mesh generation.
        P (numpy.ndarray): The initial point set (optional).
        anim (bool): Whether to animate the mesh generation.

    Returns:
        tuple: A tuple containing Node, Element, Supp, Load, and P.
    """
    if P is None:
        P = PolyMshr_RndPtSet(NElem, Domain)

    NElem = P.shape[0]
    Tol = 5e-6  # user defined tolerance for small edges
    It = 0
    Err = 1
    c = 1.5  # constant of proportionality used for calculation of 'Alpha' should be greater than 1
    BdBox = Domain.BdBox
    PFix = np.array(Domain.PFix).reshape((-1, 2))
    Area = (BdBox[1] - BdBox[0]) * (BdBox[3] - BdBox[2])
    Pc = P.copy()
    bar = progress.Bar(MaxIter, enabled=True)
    while It <= MaxIter and Err > Tol:
        Alpha = c * np.sqrt(
            Area / NElem
        )  # a distance value proportional to the width of an element
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
        # print(f"It: {It}   Error: {Err}")
        bar.increment(1, Err)
        It += 1

        if (anim == True) & (NElem <= 2000):
            PolyMshr_PlotMsh(Node, Element, NElem)

    Node, Element = PolyMshr_ExtrNds(Node, Element[:NElem])

    Node, Element = PolyMshr_CllpsEdgs(Node, Element, 0.1)

    Node, Element = PolyMshr_RsqsNds(Node, Element, NElem)

    BC = Domain.BndryCnds(Node)
    Supp = BC[0]
    Load = BC[1]

    bar.done()

    PolyMshr_PlotMsh(Node, Element, NElem, Supp, Load, wait=True)

    return Node, Element, Supp, Load, P


def PolyMshr_RndPtSet(NElem, Domain):
    """
    Generate an initial random point set of size 'NElem' for polygon mesh generation.

    Args:
        NElem (int): The number of points to generate.
        Domain (function): The domain in which points are generated.

    Returns:
        numpy.ndarray: A 2D array containing the generated points.
    """
    P = np.zeros((NElem, 2))
    BdBox = Domain.BdBox
    Ctr = 0

    while Ctr < NElem:
        Y = np.random.rand(NElem, 2)
        Y[:, 0] = (BdBox[1] - BdBox[0]) * Y[:, 0] + BdBox[0]
        Y[:, 1] = (BdBox[3] - BdBox[2]) * Y[:, 1] + BdBox[2]
        d = Domain.DistFnc(Y)
        I = np.where(d[:, -1] < 0)[0]
        NumAdded = min(NElem - Ctr, len(I))
        P[Ctr : Ctr + NumAdded, :] = Y[I[0:NumAdded], :]
        Ctr += NumAdded

    return P


def PolyMshr_FixedPoints(P, R_P, PFix):
    """
    Adjust points based on fixed points to maintain mesh quality.

    Args:
        P (numpy.ndarray): Original points.
        R_P (numpy.ndarray): Reflected points.
        PFix (numpy.ndarray): Fixed points.

    Returns:
        tuple: A tuple containing adjusted points P and R_P.
    """
    PP = np.vstack((P, R_P))
    for i in range(PFix.shape[0]):
        B, I = np.sort(
            np.sqrt((PP[:, 0] - PFix[i, 0]) ** 2 + (PP[:, 1] - PFix[i, 1]) ** 2)
        ), np.argsort(
            np.sqrt((PP[:, 0] - PFix[i, 0]) ** 2 + (PP[:, 1] - PFix[i, 1]) ** 2)
        )
        for j in range(1, 4):
            n = PP[I[j], :] - PFix[i, :]
            n = n / np.linalg.norm(n)
            PP[I[j], :] = PP[I[j], :] - n * (B[j] - B[0])

    P = PP[: P.shape[0], :]
    R_P = PP[P.shape[0] :, :]
    return P, R_P


def PolyMshr_Rflct(P, NElem, Domain, Alpha):
    """
    Reflect points at the boundary for mesh generation.

    Args:
        P (numpy.ndarray): Original points.
        NElem (int): Number of elements.
        Domain (function): The domain for boundary checks.
        Alpha (float): The propotional distance value.

    Returns:
        numpy.ndarray: Reflected points.
    """
    eps = 1e-8  # Small positive number for numerical differentiation
    # A specified parameter (0<eta<1) to adjust for numerical errors (round-off and numerical differentiation)
    eta = 0.9
    d = Domain.DistFnc(P)
    NBdrySegs = d.shape[1] - 1

    # The gradient of the distance function is computed by means of numerical differentiation
    n1 = ((Domain.DistFnc(P + np.array([[eps, 0]] * NElem))) - d) / eps
    n2 = ((Domain.DistFnc(P + np.array([[0, eps]] * NElem))) - d) / eps
    I = np.abs(d[:, 0:NBdrySegs]) < Alpha
    P1 = np.broadcast_to(P[:, 0][:, np.newaxis], (P[:, 0].shape[0], NBdrySegs))
    P2 = np.broadcast_to(P[:, 1][:, np.newaxis], (P[:, 1].shape[0], NBdrySegs))

    P1 = P1[I] - 2 * n1[:, 0:NBdrySegs][I] * d[:, 0:NBdrySegs][I]
    P2 = P2[I] - 2 * n2[:, 0:NBdrySegs][I] * d[:, 0:NBdrySegs][I]
    R_P = np.vstack((P1, P2)).T
    d_R_P = Domain.DistFnc(R_P)
    J = (np.abs(d_R_P[:, -1]) >= eta * np.abs(d[:, 0:NBdrySegs][I])) & (
        d_R_P[:, -1] > 0
    )

    R_P = np.unique(R_P[J, :], axis=0)
    return R_P


def PolyMshr_CntrdPly(Element, Node, NElem):
    """
    Compute centroids and areas for elements in the mesh.

    Args:
        Element (list): List of element vertices.
        Node (numpy.ndarray): Node coordinates.
        NElem (int): Number of elements to consider.

    Returns:
        tuple: A tuple containing centroids and areas.
    """
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

            # Compute the area of the polygon using the shoelace formula
            vx = polygon_vertices[:, 0]
            vy = polygon_vertices[:, 1]
            vxs = np.roll(vx, 1)
            vys = np.roll(vy, 1)
            temp = vx * vys - vy * vxs
            area = 0.5 * np.sum(temp)
            areas.append(area)

            # Compute the centroid of the polygon
            centroid = (
                1
                / (6 * area)
                * np.array([np.sum((vx + vxs) * temp), np.sum((vy + vys) * temp)])
            )
            centroids.append(centroid)

            counter += 1

    return np.array(centroids), np.array(areas)


def PolyMshr_ExtrNds(Node0, Element0):
    """
    Extract unique nodes and rebuild node and element lists.

    Args:
        Node0 (numpy.ndarray): Original node coordinates.
        Element0 (list): List of element vertices.

    Returns:
        tuple: A tuple containing updated Node and Element.
    """
    unique_nodes = np.unique(np.concatenate(Element0))
    cNode = np.arange(len(Node0))
    cNode[~np.in1d(cNode, unique_nodes)] = np.max(unique_nodes)
    Node, Element = PolyMshr_RbldLists(Node0, Element0, cNode)
    return Node, Element


def PolyMshr_CllpsEdgs(Node0, Element0, Tol):
    """
    Collapse small edges based on a specified tolerance.

    Args:
        Node0 (numpy.ndarray): Node coordinates.
        Element0 (list): List of element vertices.
        Tol (float): Tolerance for edge collapse.

    Returns:
        tuple: A tuple containing updated Node and Element.
    """
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
    """
    Rearrange nodes to improve mesh quality using (RCM) algorithm.

    Args:
        Node0 (numpy.ndarray): Original node coordinates.
        Element0 (list): List of element vertices.
        NElem (int): Number of elements to consider.

    Returns:
        tuple: A tuple containing updated Node and Element.
    """
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

    K = csr_matrix((s, (i, j)), shape=(NNode0, NNode0))
    p = csgraph.reverse_cuthill_mckee(K)

    cNode = np.arange(0, NNode0)
    cNode[p[:NNode0]] = np.arange(0, NNode0)

    Node, Element = PolyMshr_RbldLists(Node0, Element0, cNode)

    return Node, Element


def PolyMshr_RbldLists(Node0, Element0, cNode):
    """
    Rebuild node and element lists based on node mapping.

    Args:
        Node0 (numpy.ndarray): Original node coordinates.
        Element0 (list): List of element vertices.
        cNode (numpy.ndarray): Node mapping.

    Returns:
        tuple: A tuple containing updated Node and Element.
    """
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
    """
    Plot the polygon mesh with optional boundary conditions.

    Args:
        Node (numpy.ndarray): Node coordinates.
        Element (list): List of element vertices.
        NElem (int): Number of elements to consider.
        Supp (numpy.ndarray): Boundary support conditions (optional).
        Load (numpy.ndarray): Boundary loads (optional).
        wait (bool): Flag to pause for user interaction.

    Returns:
        None
    """
    plt.clf()

    Element = Element[:NElem]
    Node_set = set()
    for polygon in Element:
        vx = [Node[i, 0] for i in polygon]
        vy = [Node[i, 1] for i in polygon]
        Node_set.update(polygon)
        plt.fill(vx, vy, "-w", edgecolor="black")

    Node_set = Node[list(Node_set)]
    plt.plot(Node_set[:, 0], Node_set[:, 1], "bo", markersize=4)

    if Supp is not None and len(Supp) > 0:  # Plot Supp BC if specified
        plt.scatter(
            Node[Supp[:, 0], 0],
            Node[Supp[:, 0], 1],
            s=[100.0] * Supp.shape[0],
            marker=">",
        )

    if Load is not None and len(Load) > 0:  # Plot Load BC if specified
        plt.scatter(
            Node[Load[:, 0].astype(int), 0],
            Node[Load[:, 0].astype(int), 1],
            s=[100.0] * Load.shape[0],
            marker="^",
        )

    if not wait:
        plt.show(block=False)
        plt.pause(3e-1)
    else:
        plt.show()


def mesh_assessment(Node, Element, domain_area=0):
    """
    Assesses the quality of a mesh based on element aspect ratio and element area.

    This function calculates the following mesh quality metrics:

    * Maximum aspect ratio (AR) of all elements
    * Average aspect ratio of all elements
    * Average edge length across all elements
    * Range of element areas (minimum and maximum)
    * Standard deviation of element areas

    Args:
        Node (numpy.ndarray): Node coordinates.
        Element (list): List of element vertices.

    Returns:
        dict: A dictionary containing the calculated mesh quality metrics.
            - "Max. Mesh ar": Maximum aspect ratio of all elements.
            - "Average Mesh ar": Average aspect ratio of all elements.
            - "Avg. Length": Average edge length across all elements.
            - "Range of Areas": Tuple containing the minimum and maximum element areas.
            - "Standard Deviation of Areas": Standard deviation of element areas.

    Prints:
        The calculated mesh quality metrics to the console for immediate feedback.
    """
    assessment = {}

    mesh_ar = []
    all_lengthes = []
    for elem in Element:
        lengthes = []
        elem2 = np.roll(elem, 1)
        for e1, e2 in zip(elem, elem2):
            lengthes.append(np.sqrt(np.sum((Node[e1] - Node[e2]) ** 2)))
        mesh_ar.append(max(lengthes) / min(lengthes))
        all_lengthes += lengthes

    max_mesh_ar, mean_mesh_ar = max(mesh_ar), sum(mesh_ar) / len(mesh_ar)
    avg_length = sum(all_lengthes) / len(all_lengthes)

    assessment["Max. Mesh AR"], assessment["Average Mesh AR"] = (
        max_mesh_ar,
        mean_mesh_ar,
    )
    assessment["Avg. Length"] = avg_length

    print(f"Max. Aspect Ratio : {max_mesh_ar}")
    print(f"Avg. Aspect Ratio : {mean_mesh_ar}")
    print(f"Avg. Length : {avg_length}")

    areas = []

    for elem in Element:
        polygon_vertices = Node[elem, :]

        vx = polygon_vertices[:, 0]
        vy = polygon_vertices[:, 1]
        vxs = np.roll(vx, -1)
        vys = np.roll(vy, -1)
        temp = vx * vys - vy * vxs
        area = 0.5 * np.sum(temp)
        areas.append(area)

    range_areas = (min(areas), max(areas))
    standard_deviation = np.std(areas)
    assessment["Range of Areas"], assessment["Standard Deviation of Areas"] = (
        range_areas,
        standard_deviation,
    )

    print(f"Range of Areas : {range_areas[0]} - {range_areas[1]}")
    print(f"Standard Deviation of Areas : {standard_deviation}")

    return assessment
