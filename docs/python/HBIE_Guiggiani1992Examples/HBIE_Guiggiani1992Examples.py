from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
from LagrangePolynomial import *


def angle_between_two_vectors(a, b):
    return np.arccos(a@b / (np.linalg.norm(a) * np.linalg.norm(b)))

#
# Equation of external contour in polar coordinates
#     rho = rho(theta).
# see Figure 4 in Guiggiani.
# This version is based on law of sines.
#
def equation_of_external_contour_polar_coords(theta_0, theta, eta, v0, v1):
    r0  = v0 - eta
    A   = theta - theta_0
    B   = angle_between_two_vectors(v0 - v1, r0)
    C   = np.pi - A - B
    c   = np.linalg.norm(r0)
    k   = c / np.sin(C)
    rho = k * np.sin(B)
    return rho


def PolynomialBasis2d(xi1,xi2):
    return np.array([1,xi1,xi2,xi1*xi2])


class MeshGen41():
    def __init__(self):
        return

    def position(self,xi1,xi2):
        ypoints  = np.transpose(np.array([[-1,-1],[1.5,-1],[-1,1],[0.5,1]]))
        xipoints = np.transpose(np.array([[-1,-1],[1,-1],[-1,1],[1,1]]))
        A = np.zeros((4,4))
        for i in range(4):
            A[i,:] = PolynomialBasis2d(xipoints[0,i],xipoints[1,i])
        coeffs= np.linalg.solve(A,np.transpose(ypoints))
        eval = np.transpose(coeffs)@PolynomialBasis2d(xi1,xi2)

        pos = np.array([eval[0],0,eval[1]])
        return pos

class MeshGen42():
    def __init__(self):
        return

    def position(self,xi1,xi2):

        angle = np.pi/2*(xi2+1)/2
        x = np.cos(angle) 
        y = xi1 + 1
        z = np.sin(angle) 
        pos = np.array([x,y,z])
        return pos

    def derivatives(self,xi1,xi2):
        x_xi1    =  np.array([0,1,0])
        x_xi2    =  np.array([-np.pi*np.sin(np.pi*(xi2/4 + 1/4))/4,0,np.pi*np.cos(np.pi*(xi2/4 + 1/4))/4])
        x_xi1xi1 =  np.array([0,0,0])
        x_xi1xi2 =  np.array([0,0,0])
        x_xi2xi2 =  np.array([-np.pi**2*np.cos(np.pi*(xi2/4 + 1/4))/16,0,-np.pi**2*np.sin(np.pi*(xi2/4 + 1/4))/16])
        return x_xi1,x_xi2,x_xi1xi1,x_xi1xi2,x_xi2xi2

    def jacobian(self,xi1,xi2):
        jac     =  np.array([np.pi*np.cos(np.pi*(xi2/4 + 1/4))/4, 0, np.pi*np.sin(np.pi*(xi2/4 + 1/4))/4])
        return jac

    def jacobian_derivatives(self,xi1,xi2):
        jac_xi1 =  np.array([0, 0, 0])
        jac_xi2 =  np.array([-np.pi**2*np.sin(np.pi*(xi2/4 + 1/4))/16, 0, np.pi**2*np.cos(np.pi*(xi2/4 + 1/4))/16])
        return jac_xi1,jac_xi2

N = 3
n = 8
dim = 1
example = ["4.2","a"]
if example[0]=="4.1":
    mesh = MeshGen41()
    if example[1] =="a":
        eta     = np.array([0,0])
        Iexact  = -5.749237
    if example[1] =="b":
        eta     = np.array([0.66,0])
        Iexact  = -9.154585
    if example[1] =="c":
        eta     = np.array([2.0*0.885764071856287-1.0, 0.66])
        Iexact  = -15.3285
if example[0]=="4.2":
    mesh    = MeshGen42()
    if example[1] =="a":
        eta     = np.array([0,0])
        Iexact  = -0.343807
    if example[1] =="b":
        eta     = np.array([0.66,0])
        Iexact  = -0.497099
    if example[1] =="c":
        eta     = np.array([0.66,0.66])
        Iexact  = -0.877214

y = mesh.position(eta[0],eta[1])
print("y ", y)

xi1d      = np.linspace(-1,1,N)
xi1,xi2 = np.meshgrid(xi1d,xi1d)
xi = np.array([xi1.flatten(),xi2.flatten()])
x = np.zeros((3,9))

x_xi1     = np.zeros((9,3)) 
x_xi2     = np.zeros((9,3))
x_xi1xi1  = np.zeros((9,3))
x_xi2xi1  = np.zeros((9,3))
x_xi2xi2  = np.zeros((9,3))

jac       = np.zeros((9,3)) 
jac_xi1   = np.zeros((9,3)) 
jac_xi2   = np.zeros((9,3))

for i in range(9):
    x[:,i] = mesh.position(xi1.flatten()[i],xi2.flatten()[i])
    derivative = mesh.derivatives(xi1.flatten()[i],xi2.flatten()[i])
    x_xi1[i,:]    = derivative[0]
    x_xi2[i,:]    = derivative[1]
    x_xi1xi1[i,:] = derivative[2]
    x_xi2xi1[i,:] = derivative[3]
    x_xi2xi2[i,:] = derivative[4]

    jac[i,:]    = mesh.jacobian(xi1.flatten()[i],xi2.flatten()[i])
    jacobian_derivatives = mesh.jacobian_derivatives(xi1.flatten()[i],xi2.flatten()[i])
    jac_xi1[i,:]    = jacobian_derivatives[0]
    jac_xi2[i,:]    = jacobian_derivatives[1]

#edgeloop = [0,1,2,5,8,7,6,3]
edgeloop = [0,2,8,6]

#------------------------------------------------------------------------------
# Derivative matrices:
#------------------------------------------------------------------------------

# # First derivatives:
D1,D2 = Lagrange2DDerivativeMatrices(xi1d,1)

# # Second derivatives:
# D11 = D1@D1
# D12 = D1@D2
# D22 = D2@D2

# x_xi1    = D1@np.transpose(x) 
# x_xi2    = D2@np.transpose(x) 
# x_xi1xi1 = D11@np.transpose(x) 
# x_xi2xi1 = D12@np.transpose(x) 
# x_xi2xi2 = D22@np.transpose(x)
# jac = 0*x_xi1
# for i in range(x_xi1.shape[0]): 
#     jac[i,:] = np.cross(x_xi1[i,:],x_xi2[i,:])
# jac_xi1  = D1@jac
# jac_xi2  = D2@jac

interpToEta = Lagrange2DInterpMatrix(xi1d,xi1d,eta[0],eta[1]).flatten()

derivative = mesh.derivatives(eta[0],eta[1])
x_xi1_eta    = derivative[0]
x_xi2_eta    = derivative[1]
x_xi1xi1_eta = derivative[2]
x_xi2xi1_eta = derivative[3]
x_xi2xi2_eta = derivative[4]


# x_xi1_eta    = interpToEta@x_xi1   
# x_xi2_eta    = interpToEta@x_xi2   
# x_xi1xi1_eta = interpToEta@x_xi1xi1
# x_xi2xi1_eta = interpToEta@x_xi2xi1
# x_xi2xi2_eta = interpToEta@x_xi2xi2
# jac_eta      = interpToEta@jac
# jac_xi1_eta  = interpToEta@jac_xi1
# jac_xi2_eta  = interpToEta@jac_xi2

jac_eta     = mesh.jacobian(eta[0],eta[1])
derivative  = mesh.jacobian_derivatives(eta[0],eta[1])
jac_xi1_eta = derivative[0]
jac_xi2_eta = derivative[1] 

Na_eta      = interpToEta
Na_xi1_eta  = interpToEta@D1
Na_xi2_eta  = interpToEta@D2

gaussx,gaussw = np.polynomial.legendre.leggauss(n)

xx = []
yy = []
I0 = 0
Im1 = 0
Im2 = 0
for i in range(len(edgeloop)):

    i1 = edgeloop[i]
    i2 = edgeloop[(i+1)%len(edgeloop)]

    v1 = xi[:,i1]
    v2 = xi[:,i2]

    r1 = v1 - eta
    r2 = v2 - eta

    theta1  = np.arctan2(r1[1],r1[0])
    theta2  = np.arctan2(r2[1],r2[0])
    if (theta2<theta1):
        theta2 += 2*np.pi

    for theta_gx,theta_gw in zip(gaussx,gaussw):

        s     = 0.5*(theta_gx+1)
        theta = (1-s)*theta1 + s*theta2

        rho_hat = equation_of_external_contour_polar_coords(theta1,theta,eta,v1,v2)

        Ai = x_xi1_eta*np.cos(theta) + x_xi2_eta*np.sin(theta)
        Bi = 1/2*x_xi1xi1_eta*np.cos(theta)**2 + x_xi2xi1_eta*np.cos(theta)*np.sin(theta) + 1/2*x_xi2xi2_eta*np.sin(theta)**2
        Ai = Ai.flatten()
        Bi = Bi.flatten()
        A = np.linalg.norm(Ai)
        B = np.linalg.norm(Bi)
        C = np.inner(Ai,Bi)

        Ji0 = jac_eta
        Ji1 = jac_xi1_eta*np.cos(theta) + jac_xi2_eta*np.sin(theta)        

        Na0 = Na_eta
        Na1 = Na_xi1_eta*np.cos(theta) + Na_xi2_eta*np.sin(theta)

        beta = 1/A
        gamma = -np.inner(Ai,Bi)/A**4

        gi1 = Ai/A**2*(Bi@Ji0 + Ai@Ji1)

        bi0 = -Ji0
        bi1 = 3*gi1 - Ji1

        ai0 = np.outer(bi0,Na0)
        ai1 = np.outer(bi1,Na0) + np.outer(bi0,Na1)

        Sm3 = 1/A**3
        Sm2 = - 3*np.inner(Ai,Bi)/A**5

        Fm2 = -1/(4*np.pi)*Sm3*ai0
        Fm1 = -1/(4*np.pi)*(Sm2*ai0+Sm3*ai1)
        
        #Fm2 = 1/(4*np.pi)*Ji0/A**3
        #Fm1 = 1/(4*np.pi)*(-3*C*Ji0/A**5 - 3*Ai/A**5*(np.inner(Ji0,Bi) + np.inner(Ji1,Ai)) + Ji1 / A**3)

        dtheta = (theta2-theta1)/2*theta_gw
        Im2 += - Fm2*(gamma/beta**2+1/rho_hat)*dtheta
        Im1 += Fm1*np.log(rho_hat/beta)*dtheta

        for rho_gx,rho_gw in zip(gaussx,gaussw):

            s   = 0.5*(rho_gx+1)
            rho = s*rho_hat

            xi1 = eta[0] + rho*np.cos(theta)
            xi2 = eta[1] + rho*np.sin(theta)
            xxx = mesh.position(xi1,xi2)

            rvec = xxx-y
            r = np.linalg.norm(rvec)
            ri = rvec/r
            Na = Lagrange2DInterpMatrix(xi1d,xi1d,xi1,xi2).flatten()
            Ji = mesh.jacobian(xi1,xi2)
            #Ji = Na@jac
            J  = np.linalg.norm(Ji)
            ni = Ji/J
            F = - 1/(4*np.pi*r**3)*np.outer((3*ri*np.inner(ri,ni) - ni),Na)*J*rho

            drho = rho_hat/2*rho_gw
            I0 += (F - (Fm2/rho**2 + Fm1/rho))*drho*dtheta

            xx.append(xi1)
            yy.append(xi2)

I = np.sum(I0+Im1+Im2,axis=1)

if example[0]=="4.1":
    print("I0      ", I0[dim]*4*np.pi)
    print("Im1     ", Im1[dim]*4*np.pi)
    print("Im2     ", Im2[dim]*4*np.pi)
    print("I       ", I[dim]*4*np.pi)
    print("Iexact  ", Iexact)
    print("err abs ", np.abs(Iexact- I[dim]*4*np.pi))
    print("err rel ", np.abs(Iexact- I[dim]*4*np.pi)/np.abs(Iexact))

if example[0]=="4.2":
    print("I0      ", I0[2])
    print("Im1     ", Im1[2])
    print("Im2     ", Im2[2])
    print("I       ", I)
    print("Iex     ", Iexact)
    print("err abs ", np.abs(Iexact- I[2]))
    print("err rel ", np.abs(Iexact- I[2])/np.abs(Iexact))

ax = plt.axes(projection='3d')
for i in range(len(xx)):
    yyy = mesh.position(xx[i],yy[i])
    ax.plot3D(np.array([yyy[0]]),np.array([yyy[1]]),np.array([yyy[2]]), 'go')
ax.plot3D(x[0,:], x[1,:],x[2,:], 'bo')
ax.plot3D(np.array([y[0]]),np.array([y[1]]),np.array([y[2]]), 'ro')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

plt.show()
