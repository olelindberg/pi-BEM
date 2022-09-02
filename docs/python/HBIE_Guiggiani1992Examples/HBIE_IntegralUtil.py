import numpy as np

def angle_between_two_vectors(a, b):
    return np.arccos(a@b / (np.linalg.norm(a) * np.linalg.norm(b)))

#
# Equation of external contour in polar coordinates
#     rho = rho(theta).
# see Figure 4 in Guiggiani.
# This version is based on law of sines.
#
def equation_of_external_contour_polar_coords(theta_0, theta, eta, v0, v1):
    r0  = v0 - eta
    A   = theta - theta_0
    B   = angle_between_two_vectors(v0 - v1, r0)
    C   = np.pi - A - B
    c   = np.linalg.norm(r0)
    k   = c / np.sin(C)
    rho = k * np.sin(B)
    return rho