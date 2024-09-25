import numpy as np
from pyPolyMesher import PolyMesher, Domain, mesh_assessment

# # --------------------------- Mbb Domain ----------------------------------------

# from pyPolyMesher.exampleDomains import MbbDomain

# # for structured mesh
# nelx = 5
# nely = 4
# dx = 3 / nelx
# dy = 1 / nely
# x_range = np.arange(dx / 2, 3, dx)
# y_range = np.arange(dy / 2, 1, dy)
# X, Y = np.meshgrid(x_range, y_range)
# P = np.column_stack((X.ravel(), Y.ravel()))
# Node, Element, Supp, Load, P = PolyMesher(MbbDomain, 20, 30, P)
# PolyMesher(MbbDomain, 20, 30, P)

# # for unstructured mesh with given points (can be used to compare results with MATLAB code)
# P = np.array([[2.60026522, 0.03321281],
#  [0.01350291, 0.39790721],
#  [1.06513347, 0.40730687],
#  [2.54716707, 0.76359476],
#  [1.29260857, 0.39074023],
#  [2.41308243, 0.38506563],
#  [2.22468734, 0.69673122],
#  [0.95129956, 0.09203691],
#  [0.6433727,  0.67682601],
#  [2.1696533,  0.3882989 ]])
# Node, Element, Supp, Load, P = PolyMesher(MbbDomain, 10, 50, P)

# # random points
# Node, Element, Supp, Load, P = PolyMesher(MbbDomain, 50, 100, anim=True)

# # ----------------------------- Horn Domain ------------------------------------------

# from pyPolyMesher.exampleDomains import HornDomain
# HornDomain.Plot()
# Node, Element, Supp, Load, P = PolyMesher(HornDomain, 150, 50, anim=True)

# # --------------------------Wrench Domain------------------------------------------

# from pyPolyMesher.exampleDomains import WrenchDomain
# WrenchDomain.Plot()
# Node, Element, Supp, Load, P = PolyMesher(WrenchDomain, 150, 100)

# # ------------------------- Michell Domain ---------------------------------------

# from pyPolyMesher.exampleDomains import MichellDomain
# MichellDomain.Plot()
# Node, Element, Supp, Load, P = PolyMesher(MichellDomain, 20, 100)

# # with fixed points
# MichellDomain.PFix = [5, 0]
# Node, Element, Supp, Load, P = PolyMesher(MichellDomain, 20, 100)

# # ------------------------ Suspension Domain -----------------------------------------

# from pyPolyMesher.exampleDomains import SuspensionDomain
# SuspensionDomain.Plot()
# Node, Element, Supp, Load, P = PolyMesher(SuspensionDomain, 750, 150)

# # with fixed points
# SuspensionDomain.PFix = [[2, 2], [2, 16], [20, 2.5]]
# Node, Element, Supp, Load, P = PolyMesher(SuspensionDomain, 750, 150)

# ## --------------------------- Cook's Membrane Domain -----------------------------------------

# from pyPolyMesher.exampleDomains import CookDomain
# CookDomain.Plot()
# Node, Element, Supp, Load, P = PolyMesher(CookDomain, 50, 200)


### --------------------------- DXF Polygon Domain -----------------------------------------

# Import necessary functions and modules
from pyPolyMesher.dxfImporter import dxf_polygon
from pyPolyMesher.dFunctions import dPolygon

# Specify the path to the DXF file containing the polygon
dxf_file_path = 'examples/polygon1.dxf'

# Import the polygon from the DXF file
v = dxf_polygon(dxf_file_path)

# Define a Signed Distance Function (SDF) using the imported polygon
SDF = lambda P: dPolygon(P, v)

# Create a Domain object representing the DXF Polygon Domain
dxfDomain = Domain("DXF Polygon Domain", [0,100,0,100], SDF)

# Plot the DXF Polygon Domain
dxfDomain.Plot()

# Generate mesh for the DXF Polygon Domain using PolyMesher
# Specify the desired number of elements (NElem) and maximum iterations (MaxIter)
Node, Element, Supp, Load, P = PolyMesher(dxfDomain, 50, 100)

area = dxfDomain.CalculateArea()
metrics = mesh_assessment(Node, Element, area, verbose = True)

# --------------------------- How to mesh a new domain -----------------------------------------

#
# 1. Define the Bounding Box
# BdBox = [x0,y0,x1,y1]
#
# 2. Create a signed distance function
# def SDF(P):
#   # some codes to compute the signed distances of points P from domain edges
#   # you can use and combine available signed distance functions from pyPolyMesher.pydFunction module
#   return distances
#
# 3. Create a BC rule based on nodes coordinates (optional)
# def BC(nodes, BdBox):
#   # some codes to define node numbers for supported and loaded nodes
#   return [supp, load]
#
# 4. Create Domain object
# from pyPolyMesher import Domain
# NewDomain = Domain("My New Domain", BdBox, SDF, BC)
#
# 5. Plot the domain
# NewDomain.Plot()
#
# 6. Generate mesh
# from pyPolyMesher import PolyMesher
# Node, Element, Supp, Load, P = PolyMesher(NewDomain, NumberofElements, MaxIterations)
