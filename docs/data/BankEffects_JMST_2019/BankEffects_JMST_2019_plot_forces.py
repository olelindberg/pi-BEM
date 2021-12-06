import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from pandas.core.frame import DataFrame
import os
import glob

g = 9.80665
density = 1000

Lpp = 320
Sw = 9000
sheetname = 'Open data'

U=1.8013328398716324
h = 1
vel = np.sqrt(U*U)

FrL = vel/(g*Lpp)
Frh = vel/(g*h)
print('Froude number L')
print(FrL)
print('Froude number h')
print(Frh)
#meshIds = [1]

forceScale = 0.5*density*U*U*Sw;


swayforce = []
heaveforce = []
dirs = glob.glob("case3/*/")
dirs = sorted(dirs)
print(dirs)

for dir in dirs:

    filename = dir + "output/force.csv"
    print(filename)
    data = np.genfromtxt(filename, delimiter=",")
    swayforce.append(data[0][1])
    heaveforce.append(data[0][2])

heaveforce = np.array(heaveforce)
swayforce = np.array(swayforce)
#sinkage = heaveforce/(g*density*A)

print("Sway force")
print(swayforce/forceScale)

print("Sinkage")
#print(-1000*sinkage)
