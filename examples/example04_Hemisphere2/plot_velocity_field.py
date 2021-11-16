from numpy.core.numeric import NaN
import pyvista as pv
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

U   = 1
rho = 1000
g   = 9.81665

data = np.genfromtxt("/home/ole/dev/projects/pi-BEM/examples/example04d_hemisphere_in_box_hT_1.2/output/elevation.csv",delimiter=",")
xx = data[:,0]
yy = data[:,1]
zz = data[:,2]


data = np.genfromtxt("/home/ole/dev/projects/pi-BEM/examples/example04d_hemisphere_in_box_hT_1.2/output/velocity.csv",delimiter=",")

x = data[:,0]
y = data[:,1]
z = data[:,2]
pot = data[:,3]
u = -data[:,4]
v = -data[:,5]
w = -data[:,6]

V = np.sqrt(u**2 + v**2 + w**2)
r = np.sqrt(x**2 + y**2)

pot[r<=1.0] = NaN

u[r<=1.0] = NaN
v[r<=1.0] = NaN

V2 = u**2 + v**2 + w**2

pres = -rho*(-U*u + 0.5*V2)
elev =  pres/(rho*g)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(x, y, pot)


fig, ax = plt.subplots()
plt.quiver(x,y,u,v,V,cmap='jet')
ax.add_patch(plt.Circle([0,0],1,color='red',alpha=0.5))
plt.gca().set_aspect('equal','box')
plt.colorbar()


fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(x, y, elev)
ax.scatter(xx, yy, -zz,color='red')
plt.xlabel('x')
plt.ylabel('y')


plt.show()
