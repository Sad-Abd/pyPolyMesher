import numpy as np
from pyPolyMesher import *
from pyExampleDomains import *


# for structured mesh
# nelx = 5
# nely = 4
# dx = 3 / nelx
# dy = 1 / nely
# x_range = np.arange(dx / 2, 3, dx)
# y_range = np.arange(dy / 2, 1, dy)
# X, Y = np.meshgrid(x_range, y_range)
# P = np.column_stack((X.ravel(), Y.ravel()))
# PolyMesher(MbbDomain, 20, 30, P)

# for unstructured mesh with given points (can be used to compare results with MATLAB code)
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
# PolyMesher(MbbDomain, 10, 50, P)

# for unstructured mesh

# mbb_domain = MbbDomain("Mbb Domain")
# PolyMesher(mbb_domain, 50, 100)

# horn_domain = HornDomain("Horn Domain")
# PolyMesher(horn_domain, 150, 50, anim=True)

# wrench_domain = WrenchDomain("Wrench Domain")
# PolyMesher(wrench_domain, 150, 100)

# Example with Fixed Points
# michell_domain = MichellDomain("Michell Domain")
# PolyMesher(michell_domain, 20, 100)

# suspension_domain = SuspensionDomain("Suspension Domain")
# PolyMesher(suspension_domain, 750, 150)

cook_domain = CookDomain("Cook Domain")
Node, Element, Supp, Load, P = PolyMesher(cook_domain, 50, 20)
