from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
from LagrangePolynomial import *
from HBIE_Integral import *




def PolynomialBasis2d(xi1,xi2):
    return np.array([1,xi1,xi2,xi1*xi2])


class MeshGen41():
    def __init__(self):
        ypoints  = np.transpose(np.array([[-1,-1],[1.5,-1],[-1,1],[0.5,1]]))
        xipoints = np.transpose(np.array([[-1,-1],[1,-1],[-1,1],[1,1]]))
        A = np.zeros((4,4))
        for i in range(4):
            A[i,:] = PolynomialBasis2d(xipoints[0,i],xipoints[1,i])
        coeffs= np.linalg.solve(A,np.transpose(ypoints))
        self._a = coeffs[:,0]
        self._b = coeffs[:,1]
        return

    def position(self,xi1,xi2):
        x = np.transpose(self._a)@PolynomialBasis2d(xi1,xi2)
        y = np.transpose(self._b)@PolynomialBasis2d(xi1,xi2)
        pos = np.array([x,y,0])
        return pos

    def derivatives(self,xi1,xi2):
        a = self._a
        b = self._b
        x_xi1    =  np.array([a[1] + a[3]*xi2, b[1] + b[3]*xi2, 0])
        x_xi2    =  np.array([a[2] + a[3]*xi1, b[2] + b[3]*xi1, 0])
        x_xi1xi1 =  np.array([0, 0, 0])
        x_xi1xi2 =  np.array([a[3], b[3], 0])
        x_xi2xi2 =  np.array([0, 0, 0])
        return x_xi1,x_xi2,x_xi1xi1,x_xi1xi2,x_xi2xi2

    def jacobian(self,xi1,xi2):
        a = self._a
        b = self._b
        jac     =  np.array([0, 0, (a[1] + a[3]*xi2)*(b[2] + b[3]*xi1) - (a[2] + a[3]*xi1)*(b[1] + b[3]*xi2)])
        return jac

    def jacobian_derivatives(self,xi1,xi2):
        a = self._a
        b = self._b
        jac_xi1 =  np.array([0, 0, -a[3]*(b[1] + b[3]*xi2) + b[3]*(a[1] + a[3]*xi2)])
        jac_xi2 =  np.array([0, 0, a[3]*(b[2] + b[3]*xi1) + b[3]*(-a[2] - a[3]*xi1)])    
        return jac_xi1,jac_xi2


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

class LagrangePolynomial2DShapeFunction():
    def __init__(self,order):

        self._xi                  = np.linspace(-1,1,order+1)
        self._xi1,self._xi2 = np.meshgrid(self._xi,self._xi)
        self._D1,self._D2   = LagrangePolynomial2DDerivativeMatrices(self._xi,1)
        self._N = 0
        return

    def eval(self,eta):
        self._N  = LagrangePolynomial2DInterpMatrix(self._xi,self._xi,eta[0],eta[1])
        return self._N

    def derivatives(self,eta):
        return self._N@self._D1,self._N@self._D2


class PointShapeFunction():
    def __init__(self):
        return
    def eval(self,eta):
        return 1

    def derivatives(self,eta):
        return 0,0

shapefunc = PointShapeFunction()
#shapefunc = LagrangePolynomial2DShapeFunction(4)

N = 3
n = 1
example = ["4.1","a"]
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

I = HBIE_integral(mesh,shapefunc,eta,n)
I = np.sum(I,axis=1)

if example[0]=="4.1":
    print("I       ", I*4*np.pi)
    print("Iexact  ", Iexact)
    print("err abs ", np.abs(Iexact- I[2]*4*np.pi))
    print("err rel ", np.abs(Iexact- I[2]*4*np.pi)/np.abs(Iexact))

if example[0]=="4.2":
    print("I       ", I)
    print("Iex     ", Iexact)
    print("err abs ", np.abs(Iexact- I[2]))
    print("err rel ", np.abs(Iexact- I[2])/np.abs(Iexact))

plt.show()
