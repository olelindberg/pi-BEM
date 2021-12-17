
import numpy as np
from SingleTermExpansionUnknownOrder import *
from SingleTermExpansionFirstOrder import *
from SingleTermExpansionFirstAndSecondOrder import *
from SingleTermExpansionSecondOrder import *


def ErrorEstimateAndStandardDeviation(w, h, phi):
    ng = len(h)
    [phi0, alpha, p] = SingleTermExpansionUnknownOrder(w, h, phi)
    phiFit = phi0 + alpha*h**p
    errorEstimate = abs(alpha*h**p)
    standardDeviation = np.sqrt(sum(w*(phi - phiFit)**2)/((ng-1)/ng*sum(w)))

    if (p < 0.5):
        print("p is small: " + str(p))
        [std1, err1, fit1] = SingleTermExpansionFirstOrder(w, h, phi)
        [std2, err2, fit2] = SingleTermExpansionSecondOrder(w, h, phi)
        [std12, err12, fit12] = SingleTermExpansionFirstAndSecondOrder(
            w, h, phi)

        if (std1 < std2 and std1 < std12):
            standardDeviation = std1
            errorEstimate = err1
            phiFit = fit1
            p = 1.0
        elif ((std2 < std1 and std2 < std12)):
            standardDeviation = std2
            errorEstimate = err2
            phiFit = fit2
            p = 2.0
        elif ((std12 < std1 and std12 < std2)):
            standardDeviation = std12
            errorEstimate = err12
            phiFit = fit12
            p = 1.5

    elif (p > 2.0):
        print("p is large: " + str(p))

    return errorEstimate, standardDeviation, phiFit, p
