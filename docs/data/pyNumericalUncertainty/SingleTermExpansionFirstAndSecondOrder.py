
import numpy as np


def SingleTermExpansionFirstAndSecondOrder(w, h, phi):

    ng = len(h)

    wh1 = np.sum(w*h)
    wh2 = np.sum(w*h**2)
    wh3 = np.sum(w*h**3)
    wh4 = np.sum(w*h**4)

    mat = np.array([[1,   wh1, wh2],
                    [wh1, wh2, wh3],
                    [wh2, wh3, wh4]])

    vec = np.array([np.sum(w*phi),
                    np.sum(w*phi*h),
                    np.sum(w*phi*h**2)])

    tmp = np.linalg.solve(mat, vec)

    phi0 = tmp[0]
    alpha1 = tmp[1]
    alpha2 = tmp[2]

    dphi = alpha1*h + alpha2*h**2
    fit = phi0 + dphi
    std = np.sum(w*(phi-fit)**2)/(ng-3)

    return std, dphi, fit
