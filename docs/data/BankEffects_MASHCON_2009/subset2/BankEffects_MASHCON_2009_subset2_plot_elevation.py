from numpy.core.numeric import NaN
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

current_dir = os.path.dirname(os.path.abspath(__file__))

u = 0.343

rho = 1000
g = 9.80665

stagnation_pressure = 0.5*rho*u*u

data = np.genfromtxt(current_dir+"/B/output/elevation.csv", delimiter=",")
y0 = 0.265
xx = data[:, 0]
yy = data[:, 1]-y0
zz = data[:, 2]


stb = np.where(yy<=0)
x = xx[stb]
y = yy[stb]
z = zz[stb]
id = np.argsort(x)
xstb = x[id]
ystb = y[id]
zstb = z[id]


prt = np.where(yy>=0)
x = xx[prt]
y = yy[prt]
z = zz[prt]
id = np.argsort(x)
xprt = x[id]
yprt = y[id]
zprt = z[id]

hh = zstb-zprt

stb = 0.5*rho*g*np.abs(zstb)*zstb
prt = 0.5*rho*g*np.abs(zprt)*zprt

Fyprt   = np.trapz(prt,xprt)
Fystb   = np.trapz(stb,xstb)
Fy      = Fystb-Fyprt

print(Fystb)
print(Fyprt)
print(Fy)



fig = plt.figure()
#plt.plot(xprt,0*xprt + p0,'b')
plt.plot(xprt,zprt,'r')
plt.plot(xstb,zstb,'g')
plt.plot(xstb,hh,'b')
#ax.scatter(xx, yy, zz, color='red')
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y')

fig = plt.figure()
plt.plot(xprt,prt,'r')
plt.plot(xstb,stb,'g')
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y')


plt.show()
