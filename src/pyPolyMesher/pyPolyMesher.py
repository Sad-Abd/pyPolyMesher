"""
PolyMesher Module

This module provides functions for generating polygon meshes using Lloyd's algorithm.
The generated polygon meshes are suitable for various computational geometry and finite element analysis applications.

Functions:
    - PolyMesher(Domain, NElem, MaxIter, P=None, anim=False, snap_boundary=True): Perform polygon mesh generation using Lloyd's algorithm.
    - reentrant_corner_seeds(corner, e1, e2, delta): Build fixed seeds that capture a non-convex (re-entrant) corner (Talischi et al. 2012, Fig. 15).
    - assess_nodes_outside_domain(Node, SDF, tol=1e-6): Detect and quantify mesh nodes lying outside the original domain.
    - project_nodes_to_domain(Node, SDF, node_ids=None, max_step=None): Snap nodes onto the domain boundary (zero level set).
    - PolyMshr_RmvCollinearNds(Node, Element): Remove redundant nodes that are collinear within every incident element.
    - PolyMshr_RndPtSet(NElem, Domain): Generate an initial random point set of size 'NElem' for polygon mesh generation.
    - PolyMshr_FixedPoints(P, R_P, PFix): Adjust points based on fixed points to maintain mesh quality.
    - PolyMshr_Rflct(P, NElem, Domain, Alpha): Reflect points at the boundary for mesh generation.
    - PolyMshr_CntrdPly(Element, Node, NElem, P=None): Compute centroids and areas for elements in the mesh (alignment-safe).
    - PolyMshr_ExtrNds(Node0, Element0): Extract unique nodes and rebuild node and element lists.
    - PolyMshr_CllpsEdgs(Node0, Element0, Tol): Collapse small edges based on a specified tolerance.
    - PolyMshr_RsqsNds(Node0, Element0, NElem): Rearrange nodes to improve mesh quality using (RCM) algorithm.
    - PolyMshr_RbldLists(Node0, Element0, cNode): Rebuild node and element lists based on node mapping.
    - PolyMshr_PlotMsh(Node, Element, NElem, Supp=None, Load=None, wait=False): Plot the polygon mesh with optional boundary conditions.
    - mesh_assessment(Node, Element, domain_area=0): Report element-quality metrics (aspect ratio, areas, area error).

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

import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import csgraph, csr_matrix
from scipy.spatial import Voronoi, voronoi_plot_2d

from . import progress


def _boundary_node_ids(Element, n_nodes):
    """Return the ids of nodes lying on the mesh boundary.

    A mesh edge that belongs to exactly one element is a boundary edge; its two
    endpoints are boundary nodes.
    """
    from collections import Counter

    edge_count = Counter()
    for elem in Element:
        e = [int(i) for i in elem]
        m = len(e)
        for k in range(m):
            a, b = e[k], e[(k + 1) % m]
            edge_count[(a, b) if a < b else (b, a)] += 1
    ids = set()
    for (a, b), c in edge_count.items():
        if c == 1:
            ids.add(a)
            ids.add(b)
    return np.array(sorted(ids), dtype=int)


def PolyMshr_RmvCollinearNds(Node, Element, angle_tol=1e-3):
    """
    Remove redundant nodes that are collinear within every element they belong
    to (e.g. an extra node lying on a straight boundary wall between its two
    neighbours, producing two collinear edges on one wall).

    A node is removed only if it is collinear in *all* of its incident elements,
    which keeps the mesh conforming (no hanging/T-junction nodes are created).
    True corners (interior angle far from 180 deg, e.g. pinned convex corners)
    are never removed.

    Args:
        Node (numpy.ndarray): Node coordinates, shape (N, 2).
        Element (list): Element connectivity (lists/arrays of node indices).
        angle_tol (float): Relative tolerance on the cross product used to judge
            collinearity (smaller = stricter).

    Returns:
        tuple: (Node, Element) with collinear nodes removed and lists rebuilt.
    """
    from collections import defaultdict

    Elements = [[int(i) for i in e] for e in Element]
    occ = defaultdict(list)  # node -> list of (element index, position in element)
    for ei, e in enumerate(Elements):
        for pos, nd in enumerate(e):
            occ[nd].append((ei, pos))

    removable = set()
    for nd, places in occ.items():
        collinear_everywhere = True
        for ei, pos in places:
            e = Elements[ei]
            m = len(e)
            a = Node[e[(pos - 1) % m]]
            b = Node[nd]
            c = Node[e[(pos + 1) % m]]
            u = a - b
            w = c - b
            nu = np.linalg.norm(u)
            nw = np.linalg.norm(w)
            if nu < 1e-14 or nw < 1e-14:
                continue  # coincident neighbour -> treat this occurrence as redundant
            cross = u[0] * w[1] - u[1] * w[0]
            dot = u[0] * w[0] + u[1] * w[1]
            # Collinear AND between the neighbours: cross ~ 0 and vectors opposed.
            if abs(cross) > angle_tol * nu * nw or dot > 0:
                collinear_everywhere = False
                break
        if collinear_everywhere:
            removable.add(nd)

    if not removable:
        return Node, Element

    # Drop removable nodes from every element, then compact the node list.
    kept = np.array([i for i in range(Node.shape[0]) if i not in removable])
    remap = -np.ones(Node.shape[0], dtype=int)
    remap[kept] = np.arange(kept.size)
    NewNode = Node[kept]
    NewElement = []
    for e in Elements:
        ne = [remap[nd] for nd in e if nd not in removable]
        NewElement.append(np.array(ne, dtype=int))
    return NewNode, NewElement


def project_nodes_to_domain(
    Node, SDF, node_ids=None, max_iter=10, tol=1e-12, eps=1e-8, max_step=None
):
    """
    Project mesh nodes onto the domain boundary (the zero level set of the SDF).

    Boundary vertices in a reflection-based polygonal mesh do not land exactly on
    the true boundary: some poke a little outside straight walls (signed distance
    slightly positive), while reflection gaps leave others a little inside
    (signed distance slightly negative), so the meshed shape is not exactly the
    original domain (a small, persistently negative area error). This routine
    snaps nodes onto the boundary with Newton steps along the SDF gradient:

        x <- x - d(x) * grad(d) / |grad(d)|

    The same step converges to ``d = 0`` from either side.

    Args:
        Node (numpy.ndarray): Node coordinates, shape (N, 2).
        SDF (callable): Domain signed distance function (e.g. ``Domain.SDF``);
            accepts an (N, 2) array and returns an array whose last column is
            the signed distance.
        node_ids (array-like, optional): Indices of nodes to project onto the
            boundary regardless of sign (typically the boundary nodes). When
            None, only nodes that lie OUTSIDE the domain (d > tol) are pulled in.
        max_iter (int): Maximum number of Newton iterations.
        tol (float): Convergence tolerance on |d|.
        eps (float): Step used for the numerical SDF gradient.
        max_step (float, optional): If given, cap each node's total displacement
            to this distance, so a stray deep node cannot be dragged far enough
            to tangle the mesh.

    Returns:
        numpy.ndarray: A copy of Node with the selected nodes projected.
    """
    Node = np.asarray(Node, dtype=float).copy()
    ex = np.array([eps, 0.0])
    ey = np.array([0.0, eps])

    if node_ids is None:
        targeted = None
    else:
        targeted = np.asarray(node_ids, dtype=int)
        origin = Node[targeted].copy()

    for _ in range(max_iter):
        d = np.asarray(SDF(Node))[:, -1]
        if targeted is None:
            idx = np.where(d > tol)[0]
        else:
            idx = targeted[np.abs(d[targeted]) > tol]
        if idx.size == 0:
            break

        P = Node[idx]
        # Central-difference gradient of the signed distance field.
        gx = (np.asarray(SDF(P + ex))[:, -1] - np.asarray(SDF(P - ex))[:, -1]) / (2 * eps)
        gy = (np.asarray(SDF(P + ey))[:, -1] - np.asarray(SDF(P - ey))[:, -1]) / (2 * eps)
        g = np.column_stack((gx, gy))
        gnorm = np.linalg.norm(g, axis=1, keepdims=True)
        gnorm[gnorm == 0] = 1.0

        Node[idx] = P - d[idx][:, None] * g / gnorm

    if targeted is not None and max_step is not None:
        # Cap displacement of any targeted node to avoid tangling.
        disp = Node[targeted] - origin
        dist = np.linalg.norm(disp, axis=1)
        over = dist > max_step
        if over.any():
            Node[targeted[over]] = origin[over] + disp[over] * (
                max_step / dist[over][:, None]
            )

    return Node


def reentrant_corner_seeds(corner, e1, e2, delta):
    """
    Generate fixed seeds that capture a re-entrant (non-convex) corner.

    Implements the fix described by Talischi et al. (2012), "PolyMesher: a
    general-purpose mesh generator for polygonal elements written in Matlab",
    Struct Multidisc Optim, Fig. 15 and accompanying text. Near a non-convex
    corner the boundary seeds are placed independently during Lloyd's
    iterations and their reflections interfere, so the two walls meeting at the
    corner are not approximated straight. The remedy is to fix a set of seeds
    that are equidistant from the corner and all share a single common
    reflection point ``y`` placed just outside the corner (in the notch).

    Construction: an external point ``y = corner + delta * (e1 + e2)`` is placed
    in the re-entrant notch. The fixed seeds are the reflections of ``y``
    about each incident wall and through the corner point, so that each fixed
    seed reflects back onto ``y``:

        y1 = reflection of y across wall 1 (line through corner, direction e1)
        y2 = reflection of y across wall 2 (line through corner, direction e2)
        y3 = reflection of y through the corner point (point reflection)

    All three lie inside the domain at the same distance from the corner.

    Args:
        corner (array-like): The (x, y) coordinates of the re-entrant corner.
        e1 (array-like): Direction of one boundary edge leaving the corner
            (need not be unit length; it is normalized internally).
        e2 (array-like): Direction of the other boundary edge leaving the corner.
        delta (float): Offset distance of the external point ``y`` from the
            corner. Choose roughly a fraction (e.g. 0.5-1.0) of the target
            element size so the fixed seeds sit within the boundary band.

    Returns:
        list[list[float]]: The fixed seeds ``[y1, y2, y3]`` suitable for
        appending to ``Domain.PFix``.
    """
    corner = np.asarray(corner, dtype=float)
    e1 = np.asarray(e1, dtype=float)
    e2 = np.asarray(e2, dtype=float)
    e1 = e1 / np.linalg.norm(e1)
    e2 = e2 / np.linalg.norm(e2)

    # External point y placed in the re-entrant notch.
    y = corner + delta * (e1 + e2)
    v = y - corner

    def reflect_line(direction):
        # Reflect point y about the line through `corner` with unit `direction`.
        return corner + 2.0 * np.dot(v, direction) * direction - v

    y1 = reflect_line(e1)
    y2 = reflect_line(e2)
    y3 = 2.0 * corner - y  # reflection through the corner point

    return [y1.tolist(), y2.tolist(), y3.tolist()]


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


def PolyMesher(Domain, NElem, MaxIter, P=None, anim=False, snap_boundary=True):
    """
    Perform polygon mesh generation using Lloyd's algorithm.

    Args:
        Domain (function): The domain in which the mesh is generated.
        NElem (int): The desired number of elements in the mesh.
        MaxIter (int): The maximum number of iterations for mesh generation.
        P (numpy.ndarray): The initial point set (optional).
        anim (bool): Whether to animate the mesh generation.
        snap_boundary (bool): If True (default), project any final node that lies
            slightly outside the domain back onto the boundary, removing the
            small boundary overshoots inherent to reflection-based meshing.

    Returns:
        tuple: A tuple containing Node, Element, Supp, Load, and P.
    """
    if P is None:
        P = PolyMshr_RndPtSet(NElem, Domain)

    NElem = P.shape[0]
    Tol = 5e-6
    It = 0
    Err = 1
    c = 1.5  # constant of proportionality used for calculation of 'Alpha' should be greater than 1
    BdBox = Domain.BdBox
    PFix = np.array(Domain.PFix).reshape((-1, 2))
    Area = (BdBox[1] - BdBox[0]) * (BdBox[3] - BdBox[2])
    Pc = P.copy()
    bar = progress.Bar(MaxIter, enabled=True)

    allErr = []
    Node = None
    Element = None
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

        Pc, A = PolyMshr_CntrdPly(Element, Node, NElem, P)

        # Stability guard: in a correctly reflected mesh every seed-cell centroid
        # lies inside the domain. Occasionally a seed's Voronoi cell is corrupt
        # (unbounded, or bounded by near-collinear generators with far-away
        # vertices), which would fling the seed far outside and never recover,
        # leaving nodes outside the final mesh. Detect any centroid that is
        # non-finite or outside the domain and keep that seed where it is so
        # Lloyd's iteration can recover on the next pass.
        finite_mask = np.isfinite(Pc).all(axis=1)
        dPc = np.full(NElem, -1.0)
        if finite_mask.any():
            dPc[finite_mask] = Domain.DistFnc(Pc[finite_mask])[:, -1]
        bad = (~finite_mask) | (dPc > Tol)
        if bad.any():
            Pc[bad] = P[bad]
            A[bad] = 0.0

        Area = np.sum(np.abs(A))
        Err = (
            np.sqrt(np.sum((A**2) * np.sum((Pc - P) * (Pc - P), axis=1)))
            * NElem
            / (Area**1.5)
        )
        # print(f"It: {It}   Error: {Err}")
        allErr.append(Err)
        bar.increment(1, Err)
        It += 1

        if anim and (NElem <= 2000):
            PolyMshr_PlotMsh(Node, Element, NElem, Pfix=Domain.PFix)

    if Element is not None:
        Node, Element = PolyMshr_ExtrNds(Node, Element[:NElem])

        Node, Element = PolyMshr_CllpsEdgs(Node, Element, 0.1)

        Node, Element = PolyMshr_RsqsNds(Node, Element, NElem)

        if snap_boundary:
            # Snap every boundary node exactly onto the domain boundary so the
            # meshed shape matches the original (removes both outside overshoots
            # and inside reflection-gap notches that cause the area error).
            bnd_ids = _boundary_node_ids(Element, Node.shape[0])
            elem_size = np.sqrt(
                (BdBox[1] - BdBox[0]) * (BdBox[3] - BdBox[2]) / NElem
            )
            Node = project_nodes_to_domain(
                Node, Domain.SDF, node_ids=bnd_ids, max_step=elem_size
            )

            # Snapping a reflection-gap notch vertex onto a straight wall can
            # leave it collinear with its two boundary neighbours (3 colinear
            # nodes / 2 colinear edges on one wall). Drop such redundant nodes.
            Node, Element = PolyMshr_RmvCollinearNds(Node, Element)

    BC = Domain.BndryCnds(Node)
    Supp = BC[0]
    Load = BC[1]

    bar.done()

    PolyMshr_PlotMsh(Node, Element, NElem, Supp, Load, wait=True, Pfix=Domain.PFix)

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
        B, I = (
            np.sort(
                np.sqrt((PP[:, 0] - PFix[i, 0]) ** 2 + (PP[:, 1] - PFix[i, 1]) ** 2)
            ),
            np.argsort(
                np.sqrt((PP[:, 0] - PFix[i, 0]) ** 2 + (PP[:, 1] - PFix[i, 1]) ** 2)
            ),
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


def PolyMshr_CntrdPly(Element, Node, NElem, P=None):
    """
    Compute centroids and areas for elements in the mesh.

    Args:
        Element (list): List of element vertices.
        Node (numpy.ndarray): Node coordinates.
        NElem (int): Number of elements to consider.
        P (numpy.ndarray, optional): Current seed positions, shape (NElem, 2).
            Used as a safe fallback centroid for any degenerate/unbounded seed
            cell so the seed stays in place instead of being flung outside.

    Returns:
        tuple: A tuple containing centroids and areas.
    """
    # IMPORTANT: the centroid of seed `el` must end up in row `el`. We therefore
    # iterate over exactly the first NElem regions (the seed cells), producing
    # one centroid per seed in order. The previous implementation skipped
    # degenerate regions while advancing a shared counter, which let a later
    # reflected-point region fill a seed's slot and silently MISALIGNED the
    # centroids with the seeds -- scrambling seed positions and pushing seeds
    # outside the domain (unbounded cells, far Voronoi vertices, broken mesh).
    centroids = np.zeros((NElem, 2))
    areas = np.zeros(NElem)

    for el in range(NElem):
        vertices = Element[el]
        # Drop scipy's infinite-vertex marker (-1) for unbounded cells; use only
        # the finite vertices. A well-reflected seed cell is bounded and has no -1.
        finite = [v for v in vertices if v != -1]

        if len(finite) < 3:
            # Degenerate/unbounded seed cell. Do NOT average its finite vertices:
            # an unbounded cell's finite vertices can lie arbitrarily far away,
            # which would launch the seed to infinity. Keep the seed where it is
            # (if known) so Lloyd's iteration recovers on the next pass.
            if P is not None:
                centroids[el] = P[el]
            elif len(finite) > 0:
                centroids[el] = Node[finite].mean(axis=0)
            continue

        # Extract the vertex coordinates
        polygon_vertices = Node[finite]

        # Compute the area of the polygon using the shoelace formula
        vx = polygon_vertices[:, 0]
        vy = polygon_vertices[:, 1]
        vxs = np.roll(vx, 1)
        vys = np.roll(vy, 1)
        temp = vx * vys - vy * vxs
        area = 0.5 * np.sum(temp)
        areas[el] = area

        if area != 0:
            centroids[el] = (
                1
                / (6 * area)
                * np.array([np.sum((vx + vxs) * temp), np.sum((vy + vys) * temp)])
            )
        elif P is not None:
            centroids[el] = P[el]
        elif len(finite) > 0:
            centroids[el] = polygon_vertices.mean(axis=0)

    return centroids, areas


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
            if ele is None or len(ele) < 4:
                continue  # Cannot collapse triangles
            ele_array = np.asarray(ele, dtype=int)
            vx = Node0[ele_array, 0]
            vy = Node0[ele_array, 1]
            nv = len(vx)
            beta = np.arctan2(vy - np.sum(vy) / nv, vx - np.sum(vx) / nv)
            beta = np.mod(
                beta[np.roll(np.arange(len(beta)), shift=-1)] - beta, 2 * np.pi
            )
            beta_ideal = 2 * np.pi / len(ele_array)
            Edge = np.column_stack((ele_array, np.roll(ele_array, shift=-1)))
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
    Element = [[]] * len(Element0)
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


def PolyMshr_PlotMsh(Node, Element, NElem, Supp=None, Load=None, wait=False, Pfix=[]):
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
    eps = 0.1 * np.sqrt(
        (max(Node_set[:, 0]) - min(Node_set[:, 0]))
        * (max(Node_set[:, 1]) - min(Node_set[:, 1]))
        / Node.shape[0]
    )
    col = []
    Pfix_arr = np.asarray(Pfix)
    if Pfix_arr.size > 0:
        Pfix_arr = Pfix_arr.reshape((-1, 2))
        for node in Node_set:
            # Calculate distances from current node to all fixed points
            # Reshape node to (1,2) for broadcasting
            distances = np.linalg.norm(Pfix_arr - node.reshape(1, 2), axis=1)

            # Check if any fixed point is closer than threshold
            if np.min(distances) < eps:
                col.append("red")
            else:
                col.append("blue")
    else:
        col = ["blue" for each in Node_set]
    plt.scatter(Node_set[:, 0], Node_set[:, 1], c=col, s=16)
    # plt.plot(Node_set[:, 0], Node_set[:, 1], "bo", markersize=4)

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


def assess_nodes_outside_domain(Node, SDF, tol=1e-6, verbose=True):
    """
    Detect mesh nodes that fall OUTSIDE the original domain, by evaluating the
    domain's signed distance function (SDF) at every node.

    A node is "outside" when its signed distance is positive (> tol). The size of
    that positive value is exactly how far the node lies beyond the true
    boundary, so it serves as a direct geometric error measure (same units as the
    domain coordinates).

    Args:
        Node (numpy.ndarray): Node coordinates, shape (N, 2).
        SDF (callable): The domain signed distance function, e.g.
            ``Domain.SDF`` or ``Domain.DistFnc``. It must accept an (N, 2) array
            and return an array whose last column is the signed distance.
        tol (float): Distance above which a node is counted as outside.
            Default 1e-6 (ignores boundary nodes sitting at d ~ 0).
        verbose (bool): Print a summary if True.

    Returns:
        dict: {
            "num_outside": int,            # how many nodes are outside
            "frac_outside": float,         # fraction of all nodes
            "max_outside": float,          # worst (largest) positive distance
            "mean_outside": float,         # mean positive distance over outside nodes
            "sum_outside": float,          # total positive distance (severity)
            "outside_ids": numpy.ndarray,  # indices of offending nodes
            "distances": numpy.ndarray,    # signed distance for every node
        }
    """
    Node = np.asarray(Node, dtype=float)
    d = np.asarray(SDF(Node))[:, -1]

    outside_mask = d > tol
    outside_ids = np.where(outside_mask)[0]
    d_out = d[outside_mask]

    result = {
        "num_outside": int(outside_ids.size),
        "frac_outside": float(outside_ids.size) / max(1, Node.shape[0]),
        "max_outside": float(d_out.max()) if d_out.size else 0.0,
        "mean_outside": float(d_out.mean()) if d_out.size else 0.0,
        "sum_outside": float(d_out.sum()) if d_out.size else 0.0,
        "outside_ids": outside_ids,
        "distances": d,
    }

    if verbose:
        print(
            f"Nodes outside domain: {result['num_outside']} / {Node.shape[0]} "
            f"({100 * result['frac_outside']:.2f}%)"
        )
        print(f"  max  outside distance: {result['max_outside']:.6f}")
        print(f"  mean outside distance: {result['mean_outside']:.6f}")
        print(f"  sum  outside distance: {result['sum_outside']:.6f}")
        if result["num_outside"]:
            worst = outside_ids[np.argsort(-d_out)][:5]
            print("  worst offenders (id, x, y, distance):")
            for i in worst:
                print(f"    {i:5d}  ({Node[i,0]:.4f}, {Node[i,1]:.4f})  d={d[i]:.6f}")

    return result


def mesh_assessment(Node, Element, domain_area=0, verbose=True):
    """
    Assesses the quality of a mesh based on element aspect ratio and element area.

    This function calculates the following mesh quality metrics:
    * Maximum aspect ratio (AR) of all elements
    * Average aspect ratio of all elements
    * Average edge length across all elements
    * Range of element areas (minimum and maximum)
    * Standard deviation of element areas
    * Total area error between domain area and total element areas

    Args:
    Node (numpy.ndarray): Node coordinates.
    Element (list): List of element vertices.
    domain_area (float): Area of the domain (optional).
    verbose (boolean): Print mesh quality metrics (optional).

    Returns:
    dict: A dictionary containing the calculated mesh quality metrics.
    - "Max. Mesh AR": Maximum aspect ratio of all elements.
    - "Average Mesh AR": Average aspect ratio of all elements.
    - "Avg. Length": Average edge length across all elements.
    - "Range of Areas": Tuple containing the minimum and maximum element areas.
    - "Standard Deviation of Areas": Standard deviation of element areas.
    - "Total Area Error (%)": Total area error between domain area and total element areas

    Prints:
    The calculated mesh quality metrics to the console for immediate feedback.
    """
    assessment = {}
    mesh_ar = []
    all_lengths = []
    areas = []

    for elem in Element:
        elem_nodes = Node[elem]
        elem_nodes_rolled = np.roll(elem_nodes, 1, axis=0)

        # Calculate edge lengths
        edge_vectors = elem_nodes - elem_nodes_rolled
        lengths = np.sqrt(np.sum(edge_vectors**2, axis=1))

        mesh_ar.append(np.max(lengths) / np.min(lengths))
        all_lengths.extend(lengths)

        # Calculate element area
        vx, vy = elem_nodes[:, 0], elem_nodes[:, 1]
        vxs, vys = np.roll(vx, -1), np.roll(vy, -1)
        area = 0.5 * np.abs(np.sum(vx * vys - vy * vxs))
        areas.append(area)

    assessment["Max. Mesh AR"] = np.max(mesh_ar)
    assessment["Average Mesh AR"] = np.mean(mesh_ar)
    assessment["Shortest Length"] = np.min(all_lengths)
    assessment["Avg. Length"] = np.mean(all_lengths)
    assessment["Range of Areas"] = (np.min(areas), np.max(areas))
    assessment["Standard Deviation of Areas"] = np.std(areas)

    if domain_area:
        total_area = np.sum(areas)
        area_error = 100 * (total_area - domain_area) / domain_area
        assessment["Total Area Error (%)"] = area_error

    if verbose:
        for crit, val in assessment.items():
            print(f"{crit}: {val}")

    return assessment
