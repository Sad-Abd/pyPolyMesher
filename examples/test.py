import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import BSpline

from ezdxf import readfile
doc = readfile("Untitled4.dxf")
ents = doc.modelspace().query('SPLINE')

control_points = np.array(ents[0].control_points)[:,:-1]
degree = ents[0].dxf.degree
knots = np.array(ents[0].knots)

# Create a B-spline curve
bspline = BSpline(knots, control_points, degree, extrapolate="periodic")


from geomdl import NURBS, BSpline, utilities

curve = NURBS.Curve()
curve.degree = degree
curve.ctrlpts = np.vstack((control_points,control_points[0,:]))
curve.knotvector = utilities.generate_knot_vector(curve.degree, len(curve.ctrlpts))
print(knots)
print(curve.knotvector)
curve.delta = 0.01

curve_points = np.array(curve.evalpts)

# Generate points on the B-spline curve
t = np.linspace(0, 1, 1000)
# curve_points = bspline(t)

# Plot the B-spline curve
plt.plot(control_points[:, 0], control_points[:, 1], 'ro-', label='Control Points')
plt.plot(curve_points[:, 0], curve_points[:, 1], 'b-', label='B-spline Curve')
plt.title('B-spline Curve')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.legend()
plt.show()
