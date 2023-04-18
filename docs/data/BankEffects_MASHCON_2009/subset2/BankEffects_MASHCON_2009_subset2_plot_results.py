import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from pandas.core.frame import DataFrame
import os
import glob

current_dir = os.path.dirname(os.path.abspath(__file__))

g = 9.80665
density = 1000

Lpp = 4.01
A = 1.73856987
Sw = 2*1.36418245

sheetname = 'Open data'

df = pd.read_excel(
    '/home/ole/dev/projects/pi-BEM/docs/data/BankEffects_MASHCON_2009/open_data_bankeffects/subset2/opendata_subset2.xlsx', sheet_name=sheetname)
df = df.set_index('short name test')
print(df)
testIds = list(df.columns[1:])

h = np.array(list(df.loc['water depth'][1:]))
U = np.array(list(df.loc['forward component of speed vector'][1:]))
V = np.array(list(df.loc['transversal component of speed vector'][1:]))
swayforceEFD = list(df.loc['transversal force'][1:])
zaEFD = np.array(list(df.loc['sinkage fore'][1:]))
zfEFD = np.array(list(df.loc['sinkage aft'][1:]))

vel = np.sqrt(U*U + V*V)

FrL = vel/(g*Lpp)
Frh = vel/(g*h)
print("testIds:")
print(testIds)
print('Froude number L')
print(FrL)
print('Froude number h')
print(Frh)
#meshIds = [1]

swayforce = []
heaveforce = []
for testId in testIds[7]:

    filename = current_dir + "/" + testId + "/output/hydrodynamic_force.csv"
    print(filename)
    data = np.genfromtxt(filename, delimiter=",")
    print(data)
    swayforce.append(data[4])
    heaveforce.append(data[5])

heaveforce = np.array(heaveforce)
sinkage = (heaveforce-30)/(g*density*A)

print("Sway force")
print(swayforceEFD)
print(swayforce)

print("Sinkage")
print(0.5*(zaEFD+zfEFD))
print(-1000*sinkage)
