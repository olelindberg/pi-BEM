from numpy.core.numeric import NaN
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

current_dir = os.path.dirname(os.path.abspath(__file__))

u = 0.343
rho = 1000
g = 9.80665

p0 = 0.5*rho*u*u

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

pstb = rho*g*zstb
pprt = rho*g*zprt

Fyprt   = np.trapz(-pprt,xprt)
Fystb   = np.trapz(pstb,xstb)
Fy      = Fyprt+Fystb

print(Fy)

fig = plt.figure()
plt.plot(xprt,0*xprt + p0,'b')
plt.plot(xprt,pprt,'r')
plt.plot(xstb,pstb,'g')
#ax.scatter(xx, yy, zz, color='red')
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y')

plt.show()
