import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from pandas.core.frame import DataFrame
import os
import glob
import math

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

g = 9.80665
density = 1000

Lpp = 320
Tm = 20.8
vel = 3.0830504374726018
S0 = 16743

FrL = vel/math.sqrt(g*Lpp)
print('Froude number L')
print(FrL)

forceScale = 0.5*density*vel*vel*Lpp*Tm

dirs = glob.glob("./*/")
dirs = sorted(dirs)
print(dirs)

cnt = 0
swayforceAll = []
sinkageAll = []
for dir in dirs:

    swayforce = []
    heaveforce = []
    subdirs = glob.glob(dir + "/*/")
    subdirs = sorted(subdirs)
    print(subdirs)

    for subdir in subdirs:

        filename = subdir + "output/force.csv"
        print(filename)
        data = np.genfromtxt(filename, delimiter=",")
        swayforce.append(data[0][1])
        heaveforce.append(data[0][2])

    singkage = np.array(heaveforce)/(density*g*S0)/Tm
    heaveforce = np.array(heaveforce)
    swayforce = np.array(swayforce)/forceScale

    sinkageAll.append(singkage)
    swayforceAll.append(swayforce)

    print("Sway force")
    print(swayforce)
    plt.plot(ypos[cnt]+0*swayforce, - swayforce, 'ro')
    print("Sinkage")
    # print(-1000*sinkage)
    cnt = cnt + 1

sinkageAll = np.array(sinkageAll)
swayforceAll = np.array(swayforceAll)
print(sinkageAll)
print(swayforceAll)
dim = np.shape(swayforceAll)
print(dim)


plt.figure(1)
for j in range(0, dim[1]):
    plt.plot(ypos[0:3], - swayforceAll[:, j],
             '-o', label='mesh ' + str(j))
plt.plot(ypos[0:3], Ybarehull[0:3], 'k-o', label='EFD barehull')
plt.plot(ypos[0:3], Ypropulsion[0:3],
         'k--o', label='EFD propulsion')
plt.grid(True)
plt.legend()
plt.xlabel(r'$\Delta y/B$')
plt.ylabel(r'$Y/\frac{1}{2} \rho V^2 $')

plt.figure(2)
for j in range(0, dim[1]):
    plt.plot(ypos[0:3], - sinkageAll[:, j]*100,
             '-o', label='mesh ' + str(j))
plt.plot(ypos[0:3], sinkage_propulsion[0:3]*100, 'k-o', label='EFD barehull')
plt.plot(ypos[0:3], sinkage_barehull[0:3]*100,
         'k--o', label='EFD propulsion')
plt.grid(True)
plt.legend()
plt.xlabel(r'$\Delta y/B$')
plt.ylabel(r'$z/T_m~ [\%]$')

plt.show()
