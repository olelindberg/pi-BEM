
import numpy as np


def SingleTermExpansionSecondOrder(w, h, phi):

    ng = len(h)

    mat = np.array([[1,           np.sum(w*h**2)],
                    [np.sum(w*h**2), np.sum(w*h**4)]])

    vec = np.array([np.sum(w*phi),
                    np.sum(w*phi*h**2)])

    tmp = np.linalg.solve(mat, vec)

    phi0 = tmp[0]
    alpha = tmp[1]

    dphi = alpha*h**2
    fit = phi0 + dphi
    std = np.sum(w*(phi-fit)**2)/(ng-2)

    return std, dphi, fit
