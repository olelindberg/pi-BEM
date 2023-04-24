from numpy.core.numeric import NaN
import pyvista as pv
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

current_dir = os.path.dirname(os.path.abspath(__file__))

data = np.genfromtxt(current_dir + "/C3/output/elevation.csv",delimiter=",")
xx = data[:,0]
yy = data[:,1]
zz = data[:,2]


fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(xx, yy, zz,color='red')
ax.set_zlim(np.min(zz),np.max(zz))

plt.xlabel('x')
plt.ylabel('y')

plt.show()