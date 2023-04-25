import os
current_dir = os.path.dirname(os.path.abspath(__file__))
import numpy as np
import math 
import pathlib

# Sinkage is calculated by:
# s = Fz/(rho*g*A)

gravity = 9.80665
density = 1000

waterplanearea = 0.416666666666667
Tm = 0.15625

Fz = []
for i in range(1,8):
    filename = current_dir + '/wigleyhull_h2_0' + str(i)  + '/output/hydrodynamic_force.csv'
    data = np.genfromtxt(filename,delimiter=",")
    Fz.append(float(data[5]))

sinkage = np.array(Fz)/(waterplanearea*gravity*density)

print('vertical force:')
print(Fz)
print('sinkage [m]:')
print(sinkage)
print('sinkage [%]:')
print(sinkage/Tm*100)
