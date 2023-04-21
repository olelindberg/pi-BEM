
import numpy as np


def SingleTermExpansionUnknownOrder(w, h, phi):

    pmin = 0
    pmax = 4.99
    tol = 1e-10
    iter = 0
    itermax = 1000
    while (True):

        p = 0.5*(pmin+pmax)

        mat = np.array([[np.sum(w*h**(2*p)), - np.sum(w*h**p)],
                        [- np.sum(w*h**p), np.sum(w)]])
        vec = np.array([np.sum(w*phi),
                        np.sum(w*h**p*phi)])

        tmp = 1/(np.sum(w)*np.sum(w*h**(2*p)) -
                 np.sum(w*h**p)**2)*np.matmul(mat, vec)
        phi0 = tmp[0]
        alpha = tmp[1]
        fp = 2*alpha*np.sum((alpha*h**p + phi0 - phi)*np.log(h)*h**p*w)

        if abs(fp) < tol or iter > itermax:
            break

        if fp < 0:
            pmin = p
        else:
            pmax = p

        iter = iter + 1

    return phi0, alpha, p
