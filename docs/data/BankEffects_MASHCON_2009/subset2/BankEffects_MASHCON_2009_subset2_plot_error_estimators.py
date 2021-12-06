import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from pandas.core.frame import DataFrame

showPotentialErrorEstimators = True
showVelocityErrorEstimators = True

numBins = 128

dir = '/home/ole/dev/projects/pi-BEM/docs/data/BankEffects_JMST_2019/case3/mesh1/output/'


for level in range(3):

    filename = dir + "potentialErrorEstimatorLevel" + \
        str(level) + ".csv"
    print(filename)
    potentialErrorEst = np.genfromtxt(filename, delimiter=",")

    filename = dir + "/velocityErrorEstimatorLevel" + \
        str(level) + ".csv"
    print(filename)
    velocityErrorEst = np.genfromtxt(filename, delimiter=",")

    if (showPotentialErrorEstimators):
        plt.figure(2*level)
        plt.hist(potentialErrorEst, bins=numBins)

    if (showVelocityErrorEstimators):
        plt.figure(2*level+1)
        plt.hist(velocityErrorEst, bins=numBins)

plt.show()
