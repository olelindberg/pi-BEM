from ErrorEstimateAndStandardDeviation import *
import matplotlib.pyplot as plt


"""
L. Eça, M. Hoekstra, 
"A procedure for the estimation of the numerical uncertainty of CFD calculations based on grid refinement studies", 
Journal of Computational Physics 262 (2014) 104-130
"""
def NumericalUncertainty(h, phi, showPlots):

    # h = h(1: end)
    # phi = phi(1: end)

    w = 1+0*h

    # --------------------------------------------------------------------------
    # 1) Error estimate:
    # --------------------------------------------------------------------------
    [errorEstimate, standardDeviation, phiEval,
        p] = ErrorEstimateAndStandardDeviation(w, h, phi)

    w = (1/h)/np.sum(1/h)

    [errorEstimateWgt, standardDeviationWgt, phiEvalWgt,
        pWgt] = ErrorEstimateAndStandardDeviation(w, h, phi)

    if standardDeviation < standardDeviationWgt:
        stdDev = standardDeviation
        errEst = errorEstimate
        phiFit = phiEval
        p = p
    else:
        stdDev = standardDeviationWgt
        errEst = errorEstimateWgt
        phiFit = phiEvalWgt
        p = pWgt

    # --------------------------------------------------------------------------
    # 2) Data range parameter:
    # --------------------------------------------------------------------------
    ng = len(h)
    dataRangeParam = (max(phi) - min(phi))/(ng-1)

    # --------------------------------------------------------------------------
    # 3) Safety factor:
    # --------------------------------------------------------------------------
    if 0.5 <= p and p < 2.1 and stdDev < dataRangeParam:
        Fs = 1.25
    else:
        Fs = 3.0

    # --------------------------------------------------------------------------
    # 4) Uncertainty:
    # --------------------------------------------------------------------------
    if stdDev < dataRangeParam:
        Uphi = Fs*errEst + stdDev + abs(phi-phiFit)
    else:
        Uphi = 3*stdDev/dataRangeParam*(errEst + stdDev + abs(phi-phiFit))

    # --------------------------------------------------------------------------
    # Plots:
    # --------------------------------------------------------------------------

    if (showPlots):
        plt.figure(100)
        plt.semilogx(h, phi+Uphi, 'k')
        plt.semilogx(h, phi-Uphi, 'k')
        plt.semilogx(h, phi, 'r-o')
        plt.semilogx(h, phiFit, 'g-o')
        plt.show()

    return Uphi
