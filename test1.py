import numpy as np
from pyMbbDomain import MbbDomain
from pyPolyMesher import PolyMesher

nelx = 5
nely = 4

dx = 3 / nelx
dy = 1 / nely

x_range = np.arange(dx / 2, 3, dx)
y_range = np.arange(dy / 2, 1, dy)

X, Y = np.meshgrid(x_range, y_range)
P = np.column_stack((X.ravel(), Y.ravel()))

# for structured mesh 
PolyMesher(MbbDomain, 20, 30, P)

# for unstructured mesh
# PolyMesher(MbbDomain, 10, 30)