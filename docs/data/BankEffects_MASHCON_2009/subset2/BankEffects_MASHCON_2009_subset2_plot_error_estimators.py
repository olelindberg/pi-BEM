import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from pandas.core.frame import DataFrame

numBins = 32

for level in range(8):

    filename = "H/mesh1/output/potentialErrorEstimatorLevel" + \
        str(level) + ".csv"
    print(filename)
    potentialErrorEst = np.genfromtxt(filename, delimiter=",")

    filename = "H/mesh1/output/velocityErrorEstimatorLevel" + \
        str(level) + ".csv"
    print(filename)
    velocityErrorEst = np.genfromtxt(filename, delimiter=",")

    numVelAssigned = np.sum(velocityErrorEst > 1.0)

    print(numVelAssigned)
    if (True):
        plt.figure(2*level)
        plt.hist(potentialErrorEst, bins=numBins)

    plt.figure(2*level+1)
    plt.hist(velocityErrorEst, bins=numBins)

plt.show()
