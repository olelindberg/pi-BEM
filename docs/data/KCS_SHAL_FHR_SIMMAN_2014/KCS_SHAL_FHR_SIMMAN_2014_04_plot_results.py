import os
current_dir = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.insert(0, current_dir + '/../pyNumericalUncertainty')
from NumericalUncertainty import NumericalUncertainty
import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import argsort
import pandas as pd

saveFigures = True

showNumericalUncertainty = False

showHeaveForce = False
showPitchMoment = False
showSinkage = True
showPitch = False

numMeshes = 7
meshes = [0,1, 2, 3, 4, 5, 6]
g = 9.80665
density = 1000
scale = 52.667
Lpp_full = 230
Lpp_model = 1/scale*Lpp_full
tank_B_modelscale = 7
tank_B_fullscale = scale*tank_B_modelscale
CB = 0.651
zB = 5.91163
zG = 11.4
Sxx = 21035630.3
GM = 0.6
# Syy = 1.96861e+07
A = 6227.87
V0 = 52030

Tm_model = 0.2050
Tm_full = 10.8
df = pd.read_csv(current_dir + '/pmm_app_shal_fhr_App01_pure_surge.csv')
df_fs = pd.read_csv(current_dir + '/pmm_app_shal_fhr_App01_pure_surge_full_scale.csv')

U_fullscale = df_fs['surge velocity [m/s]'].values
h_fullscale = df_fs['depth [m]'].values
Frh_fullscale = U_fullscale/np.sqrt(g*h_fullscale)
FrL_fullscale = U_fullscale/np.sqrt(g*Lpp_full)

testNames = df['testId'].values
U = df['surge velocity [m/s]'].values
h = df['depth [m]'].values
sink_m = df['sink_m [mm]'].values
sink_a = df['sink_a [mm]'].values
sink_f = df['sink_f [mm]'].values
trim = df['trim [mm/m]'].values


Frh = U/np.sqrt(g*h)
FrL = U/np.sqrt(g*Lpp_model)

heaveforceAll = []
pitchmomentAll = []
elementSizeAll = []
numElementsAll = []
for testId in range(0, len(testNames)):

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

        filename = pathName + "/output/result_scalar_results_0.vtu"
        print(filename)
        mesh = pv.read(filename)

        mesh = mesh.compute_cell_sizes(length=False, volume=False)
        heaveforce.append(data[5])
        pitchmoment.append(data[7])
        numElements.append(len(mesh.get_array('Area')))
        elementSize.append(np.sqrt(np.sqrt(np.mean(mesh.get_array('Area')**2))))

    heaveforceAll.append(heaveforce)
    pitchmomentAll.append(pitchmoment)
    elementSizeAll.append(elementSize)
    numElementsAll.append(numElements)
print(np.array(numElementsAll))
print(np.array(elementSizeAll))

heaveforceAll = np.array(heaveforceAll)
pitchmomentAll = np.array(pitchmomentAll)
elementSizeAll = np.array(elementSizeAll)

metaCenter = zB+Sxx/V0
GML = metaCenter-zG

sinkage = heaveforceAll/(g*density*A)
pitch = pitchmomentAll/((g*density*V0*GML))

perm = argsort(Frh)

sink_unc = 0*sinkage
for i in range(sinkage.shape[0]):
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
    plt.xlim(0.1, 0.45)
    plt.ylim(0.0, 0.035)
    if saveFigures:
        plt.savefig("KCS_SHAL_FHR_SIMMAN_2014.png")

if (showPitchMoment):
    plt.figure(3)
    plt.plot(Frh[perm], pitchmomentAll[perm, cnt], 'o-', label="BEM")
    plt.xlabel(r"$Fr_h [ ]$")
    plt.ylabel(r"$M_y [Nm]$")
    plt.grid(True)
    plt.legend()

if (showPitch):
    plt.figure(4)
    plt.plot(Frh[perm], trim[perm], 'ko-', label="EFD")
    plt.plot(Frh[perm], -1000*pitch[perm, cnt], 'o-', label="BEM")
    plt.xlabel(r"$Fr_h [ ]$")
    plt.ylabel(r"$pitch [mm/m]$")
    plt.grid(True)
    plt.legend()

cnt = cnt + 1

plt.show()
