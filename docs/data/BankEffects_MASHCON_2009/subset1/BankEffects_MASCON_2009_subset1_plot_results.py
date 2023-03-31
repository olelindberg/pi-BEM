import numpy as np
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import argsort
from numpy.core.numeric import NaN
import pandas as pd
import os


current_dir = os.path.dirname(os.path.abspath(__file__))

g = 9.80665
density = 1000

Lpp = 4.01
h = 0.243
A = 1.73856987

y = [-0.92, 0, 0.2765, 0.49]
speed = np.array([0.458, 0.572, 0.687, 0.801])
FrL = speed/np.sqrt(g*Lpp)
Frh = speed/np.sqrt(g*h)

print(FrL)
print(Frh)


zf = np.array([[2.4,	 2.9	, 3.1, 3.3],
               [4.5,	 4.8	, 5.4, 5.6],
               [6.9,	 7.7	, 8.5, 9.6],
               [9.6,	 10.3,	10.4, 	NaN]])

za = np.array([[2.2,	 2.7,	2.8,	 3.1],
               [4.0,	 4.4,	5.0,	 5.2],
               [6.4,	 7.1,	7.9,	 9.2],
               [10.8, 14.0,	16.1,	NaN]])


#for meshId in range(1, 3):
for speed_id in range(1, 5):
    i = speed_id-1
    U = speed[i]
    heaveforce = []
    for position_id in range(1, 5):

        filename = current_dir + \
            "/position" + str(position_id) + \
            "speed" + str(speed_id) + "/output/hydrodynamic_force.csv"
        print(filename)
        data = np.genfromtxt(filename, delimiter=",")
        print(data)
        heaveforce.append(data[5])

    sinkage = np.array(heaveforce)/(g*density*A)

    plt.figure(1)
    plt.plot(y, sinkage*1000, 'o-', label="BEM")
    plt.plot(y, -za[i, :], 'ko-', label="EFD za at " + str(U) + " m/s")
    plt.plot(y, -zf[i, :], 'ko--', label="EFD zf at " + str(U) + " m/s")
    plt.xlabel(r"$y [m]$")
    plt.ylabel(r"$sinkage [mm]$")
    plt.grid(True)
    plt.legend()


plt.show()
