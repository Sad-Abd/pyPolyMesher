# import ezdxf.entities as entities

# def extract_lines_from_dxf(file_path):
#     lines = []

#     # Load the DXF file
#     doc = ezdxf.readfile(file_path)

#     # Access the model space of the DXF file
#     msp = doc.modelspace()
#     ents = doc.modelspace().query('CIRCLE LINE ARC POLYLINE ELLIPSE SPLINE SHAPE LWPOLYLINE')
#     for e in ents:
#         print(e)
#     # Iterate through entities in the model space
#     for entity in msp.query('LINE'):
#         start_point = (entity.dxf.start.x, entity.dxf.start.y)
#         end_point = (entity.dxf.end.x, entity.dxf.end.y)
#         lines.append((start_point, end_point))

#     return lines

# Specify the path to your DXF file


# # Get the list of lines from the DXF file
# lines = extract_lines_from_dxf(dxf_file_path)
# # Print the list of lines
# for line in lines:
#     print(f"Start Point: {line[0]}, End Point: {line[1]}")

from ezdxf import readfile
doc = readfile("Untitled4.dxf")
ents = doc.modelspace().query('SPLINE')

import numpy as np
from scipy.interpolate import BSpline

control_points = np.array(ents[0].control_points)[:,:-1]
degree = ents[0].dxf.degree
knot_vector = np.array(ents[0].knots)
domain_size = 10


def create_sdf(control_points, degree, knot_vector, domain_size):
    """
    Create a signed distance function (SDF) for a closed B-spline curve.

    Parameters:
        control_points (numpy.ndarray): Control points of the B-spline curve.
        degree (int): Degree of the B-spline polynomial.
        knot_vector (numpy.ndarray): Knot vector of the B-spline curve.
        domain_size (int): Size of the domain for SDF calculation.

    Returns:
        numpy.ndarray: SDF grid for the specified domain.
    """
    # Create a B-spline curve
    bspline = BSpline(knot_vector, control_points, degree, extrapolate="periodic")

    # Evaluate the B-spline curve at points in the domain
    points = np.linspace(0, 1, domain_size)
    curve_points = bspline(points)
    print(curve_points)
    # Calculate signed distances for each point in the domain
    sdf = np.zeros((domain_size, domain_size))
    for i in range(domain_size):
        for j in range(domain_size):
            point = np.array([points[i], points[j]])
            distances = np.linalg.norm(curve_points - point, axis=1)
            sdf[i, j] = np.min(distances)

    return sdf

print(create_sdf(control_points, degree, knot_vector, domain_size))
