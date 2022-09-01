from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
from LagrangePolynomial import *

class Sphere():
    def __init__(self):
        return

    def position(self,xi1,xi2):

        phi     = np.pi/4*(xi1+1)/2
        theta   = np.pi/2-np.pi/4*(xi2+1)/2

        x = np.cos(phi)*np.sin(theta) 
        y = np.sin(phi)*np.sin(theta)
        z = np.cos(phi)*np.cos(theta) 
        pos = np.array([x,y,z])
        return pos

    def derivatives(self,xi1,xi2):
        x_xi1    =  [-np.pi*np.sin(np.pi*(xi1/8 + 1/8))*np.cos(np.pi*(xi2/8 + 1/8))/8, np.pi*np.cos(np.pi*(xi1/8 + 1/8))*np.cos(np.pi*(xi2/8 + 1/8))/8, -np.pi*np.sin(np.pi*(xi1/8 + 1/8))*np.sin(np.pi*(xi2/8 + 1/8))/8]
        x_xi2    =  [-np.pi*np.sin(np.pi*(xi2/8 + 1/8))*np.cos(np.pi*(xi1/8 + 1/8))/8, -np.pi*np.sin(np.pi*(xi1/8 + 1/8))*np.sin(np.pi*(xi2/8 + 1/8))/8, np.pi*np.cos(np.pi*(xi1/8 + 1/8))*np.cos(np.pi*(xi2/8 + 1/8))/8]
        x_xi1xi1 =  [-np.pi**2*np.cos(np.pi*(xi1/8 + 1/8))*np.cos(np.pi*(xi2/8 + 1/8))/64, -np.pi**2*np.sin(np.pi*(xi1/8 + 1/8))*np.cos(np.pi*(xi2/8 + 1/8))/64, -np.pi**2*np.sin(np.pi*(xi2/8 + 1/8))*np.cos(np.pi*(xi1/8 + 1/8))/64]
        x_xi1xi2 =  [np.pi**2*np.sin(np.pi*(xi1/8 + 1/8))*np.sin(np.pi*(xi2/8 + 1/8))/64, -np.pi**2*np.sin(np.pi*(xi2/8 + 1/8))*np.cos(np.pi*(xi1/8 + 1/8))/64, -np.pi**2*np.sin(np.pi*(xi1/8 + 1/8))*np.cos(np.pi*(xi2/8 + 1/8))/64]
        x_xi2xi2 =  [-np.pi**2*np.cos(np.pi*(xi1/8 + 1/8))*np.cos(np.pi*(xi2/8 + 1/8))/64, -np.pi**2*np.sin(np.pi*(xi1/8 + 1/8))*np.cos(np.pi*(xi2/8 + 1/8))/64, -np.pi**2*np.sin(np.pi*(xi2/8 + 1/8))*np.cos(np.pi*(xi1/8 + 1/8))/64]
        return x_xi1,x_xi2,x_xi1xi1,x_xi1xi2,x_xi2xi2

    def jacobian(self,xi1,xi2):
        jac     =  [-np.pi**2*np.sin(np.pi*(xi1/8 + 1/8))**2*np.sin(np.pi*(xi2/8 + 1/8))**2/64 + np.pi**2*np.cos(np.pi*(xi1/8 + 1/8))**2*np.cos(np.pi*(xi2/8 + 1/8))**2/64, np.pi**2*np.sin(np.pi*(xi1/8 + 1/8))*np.sin(np.pi*(xi2/8 + 1/8))**2*np.cos(np.pi*(xi1/8 + 1/8))/64 + np.pi**2*np.sin(np.pi*(xi1/8 + 1/8))*np.cos(np.pi*(xi1/8 + 1/8))*np.cos(np.pi*(xi2/8 + 1/8))**2/64, np.pi**2*np.sin(np.pi*(xi1/8 + 1/8))**2*np.sin(np.pi*(xi2/8 + 1/8))*np.cos(np.pi*(xi2/8 + 1/8))/64 + np.pi**2*np.sin(np.pi*(xi2/8 + 1/8))*np.cos(np.pi*(xi1/8 + 1/8))**2*np.cos(np.pi*(xi2/8 + 1/8))/64]
        return jac

    def jacobian_derivatives(self,xi1,xi2):
        jac_xi1 =  [-np.pi**3*np.sin(np.pi*(xi1/8 + 1/8))*np.sin(np.pi*(xi2/8 + 1/8))**2*np.cos(np.pi*(xi1/8 + 1/8))/256 - np.pi**3*np.sin(np.pi*(xi1/8 + 1/8))*np.cos(np.pi*(xi1/8 + 1/8))*np.cos(np.pi*(xi2/8 + 1/8))**2/256, -np.pi**3*np.sin(np.pi*(xi1/8 + 1/8))**2*np.sin(np.pi*(xi2/8 + 1/8))**2/512 - np.pi**3*np.sin(np.pi*(xi1/8 + 1/8))**2*np.cos(np.pi*(xi2/8 + 1/8))**2/512 + np.pi**3*np.sin(np.pi*(xi2/8 + 1/8))**2*np.cos(np.pi*(xi1/8 + 1/8))**2/512 + np.pi**3*np.cos(np.pi*(xi1/8 + 1/8))**2*np.cos(np.pi*(xi2/8 + 1/8))**2/512, 0]
        jac_xi2 =  [-np.pi**3*np.sin(np.pi*(xi1/8 + 1/8))**2*np.sin(np.pi*(xi2/8 + 1/8))*np.cos(np.pi*(xi2/8 + 1/8))/256 - np.pi**3*np.sin(np.pi*(xi2/8 + 1/8))*np.cos(np.pi*(xi1/8 + 1/8))**2*np.cos(np.pi*(xi2/8 + 1/8))/256, 0, -np.pi**3*np.sin(np.pi*(xi1/8 + 1/8))**2*np.sin(np.pi*(xi2/8 + 1/8))**2/512 + np.pi**3*np.sin(np.pi*(xi1/8 + 1/8))**2*np.cos(np.pi*(xi2/8 + 1/8))**2/512 - np.pi**3*np.sin(np.pi*(xi2/8 + 1/8))**2*np.cos(np.pi*(xi1/8 + 1/8))**2/512 + np.pi**3*np.cos(np.pi*(xi1/8 + 1/8))**2*np.cos(np.pi*(xi2/8 + 1/8))**2/512]
        return jac_xi1,jac_xi2

N = 3
mesh = Sphere()

xi1d      = np.linspace(-1,1,N)

tmp1,tmp2 = np.meshgrid(xi1d,xi1d)
xi = np.transpose(np.array([tmp1.flatten(),tmp2.flatten()]))

x         = np.zeros((N*N,3))
x_xi1     = np.zeros((N*N,3)) 
x_xi2     = np.zeros((N*N,3))
x_xi1xi1  = np.zeros((N*N,3))
x_xi1xi2  = np.zeros((N*N,3))
x_xi2xi2  = np.zeros((N*N,3))
jac       = np.zeros((N*N,3)) 
jac_xi1   = np.zeros((N*N,3)) 
jac_xi2   = np.zeros((N*N,3))
for i in range(N*N):
    x[i,:]                  = mesh.position(xi[i,0],xi[i,1])
    derivative              = mesh.derivatives(xi[i,0],xi[i,1])
    x_xi1[i,:]              = derivative[0]
    x_xi2[i,:]              = derivative[1]
    x_xi1xi1[i,:]           = derivative[2]
    x_xi1xi2[i,:]           = derivative[3]
    x_xi2xi2[i,:]           = derivative[4]
    jac[i,:]                = mesh.jacobian(xi[i,0],xi[i,1])
    jacobian_derivatives    = mesh.jacobian_derivatives(xi[i,0],xi[i,1])
    jac_xi1[i,:]            = jacobian_derivatives[0]
    jac_xi2[i,:]            = jacobian_derivatives[1]

D1,D2 = Lagrange2DDerivativeMatrices(xi1d,1)

# Second derivatives:
D11 = D1@D1
D12 = D1@D2
D22 = D2@D2
x_xi1_h    = D1@x 
x_xi2_h    = D2@x 
x_xi1xi1_h = D11@x 
x_xi1xi2_h = D12@x 
x_xi2xi2_h = D22@x 
jac_h = 0*x_xi1_h
for i in range(N*N): 
    jac_h[i,:] = np.cross(x_xi1_h[i,:],x_xi2_h[i,:])
jac_xi1_h  = D1@jac_h
jac_xi2_h  = D2@jac_h

plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(x[:,0], x[:,1],x[:,2], 'bo')
ax.set_xlabel("x1")
ax.set_ylabel("x2")
ax.set_zlabel("x3")


for dim in range(3):
    plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(xi[:,0], xi[:,1],x[:,dim], 'bo')
    ax.set_xlabel("xi1")
    ax.set_ylabel("xi2")
    ax.set_zlabel("x"+str(dim))

    plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(xi[:,0], xi[:,1],x_xi1[:,dim], 'bo')
    ax.plot3D(xi[:,0], xi[:,1],x_xi1_h[:,dim], 'rx')
    ax.set_xlabel("xi1")
    ax.set_ylabel("xi2")
    ax.set_zlabel("x_xi1 "+str(dim))

    plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(xi[:,0], xi[:,1],x_xi2[:,dim], 'bo')
    ax.plot3D(xi[:,0], xi[:,1],x_xi2_h[:,dim], 'rx')
    ax.set_xlabel("xi1")
    ax.set_ylabel("xi2")
    ax.set_zlabel("x_xi2 "+str(dim))

    plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(xi[:,0], xi[:,1],x_xi1xi1[:,dim], 'bo')
    ax.plot3D(xi[:,0], xi[:,1],x_xi1xi1_h[:,dim], 'rx')
    ax.set_xlabel("xi1")
    ax.set_ylabel("xi2")
    ax.set_zlabel("x_xi1xi1 "+str(dim))

    plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(xi[:,0], xi[:,1],x_xi1xi2[:,dim], 'bo')
    ax.plot3D(xi[:,0], xi[:,1],x_xi1xi2_h[:,dim], 'rx')
    ax.set_xlabel("xi1")
    ax.set_ylabel("xi2")
    ax.set_zlabel("x_xi1xi2 "+str(dim))

    plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(xi[:,0], xi[:,1],x_xi2xi2[:,dim], 'bo')
    ax.plot3D(xi[:,0], xi[:,1],x_xi2xi2_h[:,dim], 'rx')
    ax.set_xlabel("xi1")
    ax.set_ylabel("xi2")
    ax.set_zlabel("x_xi2xi2 "+str(dim))


    plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(xi[:,0], xi[:,1],jac[:,dim], 'bo')
    ax.plot3D(xi[:,0], xi[:,1],jac_h[:,dim], 'rx')
    ax.set_xlabel("xi1")
    ax.set_ylabel("xi2")
    ax.set_zlabel("jac "+str(dim))

    plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(xi[:,0], xi[:,1],jac_xi1[:,dim], 'bo')
    ax.plot3D(xi[:,0], xi[:,1],jac_xi1_h[:,dim], 'rx')
    ax.set_xlabel("xi1")
    ax.set_ylabel("xi2")
    ax.set_zlabel("jac_xi1 "+str(dim))

    plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(xi[:,0], xi[:,1],jac_xi2[:,dim], 'bo')
    ax.plot3D(xi[:,0], xi[:,1],jac_xi2_h[:,dim], 'rx')
    ax.set_xlabel("xi1")
    ax.set_ylabel("xi2")
    ax.set_zlabel("jac_xi2 "+str(dim))

plt.show()
