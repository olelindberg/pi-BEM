import numpy as np
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import argsort
import pandas as pd

numMeshes = 8

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
#Syy = 1.96861e+07
A = 1.73856987
V0 = 52030
Sw = 2*1.36418245

Tm_model = 0.2050
Tm_full = 10.8

print(zG+Tm_full)

y = [-0.92, 0, 0.2765, 0.49]


heaveforce = []
pitchmoment = []
for testId in range(1, 5):

    filename = "case" + str(testId) + "/output/force.csv"
    print(filename)
    data = np.genfromtxt(filename, delimiter=",")

    heaveforce.append(data[0][2])
    pitchmoment.append(data[0][4])

heaveforce = np.array(heaveforce)
pitchmoment = np.array(pitchmoment)

U = 1.7997864228846707
forceScale = 0.5*density*U*U*Sw
CZ = heaveforce/forceScale

speed = [0.458, 0.572, 0.687, 0.801]

for U in speed:
    forceScale = 0.5*density*U*U*Sw
    heaveforce = CZ*forceScale

    sinkage = heaveforce/(g*density*A)

    metaCenter = zB+Sxx/V0
    GML = metaCenter-zG
    pitch = pitchmoment/((g*density*V0*GML))

    plt.figure(1)
    plt.plot(y, heaveforce, 'o-', label="BEM")
    plt.xlabel(r"$Fr_h$")
    plt.ylabel(r"$F_z [N]$")
    plt.grid(True)
    plt.legend()

    plt.figure(2)
    plt.plot(y, sinkage*1000, 'o-', label="BEM")
    plt.xlabel(r"$y$")
    plt.ylabel(r"$sinkage [m]$")
    plt.grid(True)
    plt.legend()

plt.figure(3)
plt.plot(y, pitchmoment, 'o-', label="BEM")
plt.xlabel(r"$Fr_h [ ]$")
plt.ylabel(r"$M_y [Nm]$")
plt.grid(True)
plt.legend()

plt.figure(4)
plt.plot(y, -1000*pitch, 'o-', label="BEM")
plt.xlabel(r"$Fr_h [ ]$")
plt.ylabel(r"$pitch [mm/m]$")
plt.grid(True)
plt.legend()

plt.show()
