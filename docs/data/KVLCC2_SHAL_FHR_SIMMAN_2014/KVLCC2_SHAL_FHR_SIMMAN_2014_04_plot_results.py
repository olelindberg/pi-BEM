import os
current_dir = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.insert(0, current_dir + '/../pyNumericalUncertainty')
from NumericalUncertainty import NumericalUncertainty
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import argsort
import pandas as pd
import glob

meshes = [2,3,4,5,6]
showNumericalUncertainty = False
showHeaveForce = False
showPitchMoment = False
showSinkage = True
showPitch = False
saveFigures = True

g = 9.80665
density = 1000
scale = 75
Lpp_full = 320
Tm_full = 20.8
A = 16963.3
V0 = 312600
zG = 20.8
zB = 0.0
Sxx = 0.0

Lpp_model = 1/scale*Lpp_full
Tm_model = 1/scale*Tm_full


df = pd.read_csv(current_dir + '/pmm_app_shal_fhr_App01_pure_surge.csv')
df_fs = pd.read_csv(current_dir + '/pmm_app_shal_fhr_App01_pure_surge_full_scale.csv')

U_fullscale = df_fs['surge velocity [m/s]'].values
h = df_fs['depth [m]'].values
Frh = U_fullscale/np.sqrt(g*h)
FrL = U_fullscale/np.sqrt(g*Lpp_full)
print(Frh)
print(FrL)


testIds = df['testId'].values
testNames = df['testId'].values
U = df['surge velocity [m/s]'].values
h = df['depth [m]'].values
sink_m = df['sink_m [mm]'].values

Frh = U/np.sqrt(g*h)
FrL = U/np.sqrt(g*Lpp_model)
print(Frh)
print(FrL)

perm = np.argsort(Frh)

heaveforceAll = []
pitchmomentAll = []
elementSizeAll = []
numElementsAll = []

print(testNames)
for testId in range(0,len(testNames)):

    testName = testNames[testId]

    heaveforce = []
    pitchmoment = []
    elementSize = []
    numElements = []

    for meshId in meshes:

        meshName = "mesh0" + str(meshId)
        pathName = current_dir + "/" + testName + "/" + meshName

        filename = pathName + "/output/hydrodynamic_force.csv"
        print(filename)
        data = np.genfromtxt(filename, delimiter=",")

        filename = pathName + "/output/result_scalar_results.vtu"
        print(filename)
        mesh = pv.read(filename)
        mesh = mesh.compute_cell_sizes(length=False, volume=False)

        heaveforce.append(data[5])
        pitchmoment.append(data[7])
        elementSize.append(np.sqrt(np.sqrt(np.mean(mesh.get_array('Area')**2))))
        numElements.append(mesh.number_of_cells)

    heaveforceAll.append(heaveforce)
    pitchmomentAll.append(pitchmoment)
    elementSizeAll.append(elementSize)
    numElementsAll.append(numElements)

heaveforceAll = np.array(heaveforceAll)
pitchmomentAll = np.array(pitchmomentAll)
elementSizeAll = np.array(elementSizeAll)

print(numElementsAll)
print(elementSizeAll)

metaCenter = zB+Sxx/V0
GML = metaCenter-zG

sinkage = heaveforceAll/(g*density*A)
pitch = pitchmomentAll/((g*density*V0*GML))

perm = argsort(Frh)

sink_unc = 0*sinkage
for i in range(sinkage.shape[1]):
    sink_unc[i, :] = NumericalUncertainty(
        elementSizeAll[i, :], sinkage[i, :], showNumericalUncertainty)

cnt = sinkage.shape[1]-1

if (showHeaveForce):
    plt.figure(1)
    plt.plot(Frh[perm], heaveforceAll[perm, cnt], 'o-', label="BEM")
    plt.xlabel(r"$Fr_h$")
    plt.ylabel(r"$F_z [N]$")
    plt.grid(True)
    plt.legend()

if (showSinkage):
    plt.figure(2)
    plt.plot(Frh[perm], 0.001*sink_m[perm]/Tm_model, 'ko-', label="EFD")
    plt.fill_between(Frh[perm],
                     -(sinkage[perm, cnt] +
                       sink_unc[perm, cnt]) / Tm_full,
                     -(sinkage[perm, cnt] -
                       sink_unc[perm, cnt]) / Tm_full, color='lightgray', label="95% confidence band")
    plt.plot(Frh[perm], -sinkage[perm, cnt] /
             Tm_full, 'ro-', label="BEM")
    plt.xlabel(r"$Fr_h$")
    plt.ylabel(r"$sinkage/T_m$")
    plt.grid(True)
    plt.legend()
    plt.xlim(0.08, 0.24)
    plt.ylim(0.0, 0.014)
    if saveFigures:
        plt.savefig("KVLCC2_SHAL_FHR_SIMMAN_2014_sinkage.png")

if (showPitchMoment):
    plt.figure(3)
    plt.plot(Frh[perm], pitchmomentAll[cnt, perm], 'o-', label="BEM")
    plt.xlabel(r"$Fr_h [ ]$")
    plt.ylabel(r"$M_y [Nm]$")
    plt.grid(True)
    plt.legend()

if (showPitch):
    plt.figure(4)
    plt.plot(Frh[perm], trim[perm], 'ko-', label="EFD")
    plt.plot(Frh[perm], -1000*pitch[cnt, perm], 'o-', label="BEM")
    plt.xlabel(r"$Fr_h [ ]$")
    plt.ylabel(r"$pitch [mm/m]$")
    plt.grid(True)
    plt.legend()

cnt = cnt + 1

plt.show()
