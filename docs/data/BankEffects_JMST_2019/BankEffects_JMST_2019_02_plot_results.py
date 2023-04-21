import os
import sys

current_dir = os.path.dirname(os.path.realpath(__file__))
print(current_dir)
sys.path.insert(0, current_dir + '/../pyNumericalUncertainty')

from NumericalUncertainty import NumericalUncertainty
import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from pandas.core.frame import DataFrame
import glob
import math


showUncertaintyPlots = True
saveFigures = True
testNames = ['case1', 'case2', 'case3']
meshes = [1, 2, 3, 4, 5]

ywall = -2.7
beam_modelscale = 0.7733
Tm_modelscale = 0.2776
sinkage_propulsion = 1/1000 * \
    np.array([3.62, 3.558, 2.783, 2.397, 3.977])/Tm_modelscale
sinkage_barehull = 1/1000 * \
    np.array([3.356, 3.229, 2.738, 2.201, 3.714])/Tm_modelscale
ypos = (np.array([-2.1825, -2.1134, -
        1.7269, -1.7269, -1.7269])-ywall-beam_modelscale/2)/beam_modelscale

Ypropulsion = [0.1087, 0.1134, 0.04782, 0.04316, 0.02994]
Ybarehull = [0.08023, 0.06349, 0.03042, 0.02569, 0.01052]

print(sinkage_propulsion[0:3])
print(sinkage_barehull[0:3])
print(ypos[0:3])
print(Ypropulsion[0:3])
print(Ybarehull[0:3])

g = 9.80665
density = 1000

Lpp = 320
Tm = 20.8
vel = 3.0830504374726018
A = 16743

FrL = vel/math.sqrt(g*Lpp)

forceScale = 0.5*density*vel*vel*Lpp*Tm

heaveforceAll = []
swayforceAll = []
elementSizeAll = []
numElementsAll = []

for testId in range(0,len(testNames)):

    testName = testNames[testId]

    swayforce = []
    heaveforce = []
    elementSize = []
    numElements = []

    for meshId in meshes:

        meshName = "mesh" + str(meshId)
        pathName = current_dir + "/" + testName + "/" + meshName

        filename = pathName + "/output/hydrodynamic_force.csv"
        data = np.genfromtxt(filename, delimiter=",")
        print(data)
        filename = pathName + "/output/result_scalar_results.vtu"
        mesh = pv.read(filename)
        mesh = mesh.compute_cell_sizes(length=False, volume=False)

        swayforce.append(data[4])
        heaveforce.append(data[5])
        elementSize.append(np.sqrt(np.min(mesh.get_array('Area'))))
        numElements.append(mesh.number_of_cells)

    swayforceAll.append(swayforce)
    heaveforceAll.append(heaveforce)
    elementSizeAll.append(elementSize)
    numElementsAll.append(numElements)

heaveforceAll = np.array(heaveforceAll)
swayforceAll = np.array(swayforceAll)
elementSizeAll = np.array(elementSizeAll)

print(elementSizeAll)
print(numElementsAll)
sinkage = heaveforceAll/(g*density*A)

sway_unc = 0*swayforceAll
for i in range(swayforceAll.shape[0]):
    sway_unc[i, :] = NumericalUncertainty(
        elementSizeAll[i, :], swayforceAll[i, :]/forceScale, showUncertaintyPlots)

sink_unc = 0*sinkage
for i in range(0, sinkage.shape[0]):
    sink_unc[i, :] = NumericalUncertainty(
        elementSizeAll[i, :], sinkage[i, :], showUncertaintyPlots)


cnt = sinkage.shape[1]-1

plt.figure(1)
plt.plot(ypos[0:3], Ybarehull[0:3], 'k-o', label='EFD barehull')
plt.plot(ypos[0:3], Ypropulsion[0:3], 'k--o', label='EFD propulsion')
plt.fill_between(ypos[0:3],
                 -(swayforceAll[:, cnt] / forceScale +
                   sway_unc[:, cnt]),
                 -(swayforceAll[:, cnt] / forceScale -
                   sway_unc[:, cnt]), color='lightgray', label="95% confidence band")
plt.plot(ypos[0:3], - swayforceAll[:, cnt]/forceScale, 'r-o', label='BEM')
plt.grid(True)
plt.legend()
plt.xlabel(r'$\Delta y/B$')
plt.ylabel(r'$Y/(\frac{1}{2} \rho V^2 L_{pp} T_m )$')
plt.xlim([0.1, 0.8])
plt.ylim([0.02, 0.14])
if saveFigures:
    plt.savefig("BankEffects_JMST_2019_swayforce.png")

plt.figure(2)
plt.plot(ypos[0:3], sinkage_propulsion[0:3], 'k-o', label='EFD barehull')
plt.plot(ypos[0:3], sinkage_barehull[0:3], 'k--o', label='EFD propulsion')
plt.fill_between(ypos[0:3],
                 -(sinkage[:, cnt] +
                   sink_unc[:, cnt]) / Tm,
                 -(sinkage[:, cnt] -
                   sink_unc[:, cnt]) / Tm, color='lightgray', label="95% confidence band")

plt.plot(ypos[0:3], - sinkage[:, cnt]/Tm, 'r-o', label='BEM')
plt.grid(True)
plt.legend()
plt.xlabel(r'$\Delta y/B$')
plt.ylabel(r'$sinkage/T_m$')
plt.xlim([0.1, 0.8])
plt.ylim([0.0095, 0.0135])
if saveFigures:
    plt.savefig("BankEffects_JMST_2019_sinkage.png")

plt.show()
