from sympy import *
from sympy.integrals.trigonometry import trigintegrate
from sympy.plotting import plot3d
from sympy.integrals.quadrature import gauss_legendre
import math


def string_replacements(strng):
    strng = str(strng).replace("(x, y, z, t)","")
    strng = str(strng).replace("Derivative","D")
    strng = str(strng).replace("1.0*","")

    #strng = str(strng).replace("r(x, y, z)","r")
    #strng = str(strng).replace("D(r, x)","rx")
    #strng = str(strng).replace("D(r, y)","ry")
    #strng = str(strng).replace("D(r, z)","rz")
    #strng = str(strng).replace("D(r, (x, 2))","rxx")
    #strng = str(strng).replace("D(r, (y, 2))","ryy")
    #strng = str(strng).replace("D(r, (z, 2))","rzz")
    #strng = str(strng).replace("D(r, x, y)","rxy")
    #strng = str(strng).replace("D(r, x, z)","rxz")
    #strng = str(strng).replace("D(r, y, z)","ryz")
    #strng = str(strng).replace("((x - x0)**2 + (y - y0)**2 + (z - z0)**2)**(-0.5)","1.0/r")
    #strng = str(strng).replace("((x - x0)**2 + (y - y0)**2 + (z - z0)**2)**(-1.5)","1.0/r3")
    
    #strng = str(strng).replace("*1.0","")
    #strng = str(strng).replace("(x - x0)","dx")
    #strng = str(strng).replace("(y - y0)","dy")
    #strng = str(strng).replace("(z - z0)","dz")
    #strng = str(strng).replace("**2","2")
    #strng = str(strng).replace("**3","3")
    #strng = str(strng).replace("2*pi","two_pi")
    #strng = str(strng).replace("4*pi","four_pi")
    return strng

def expandEq(eq):

    eq = expand(eq)
    eq = collect(eq,e)
    eq0 = eq.coeff(e,0)
    eq1 = eq.coeff(e,1)

    return eq0,eq1


nx,ny,nz,x,y,z,x0,y0,z0,pi,t,e,U,g = symbols('nx ny nz x y z x0 y0 z0 pi t e U g')

eta0 = Function('eta0')(x,y,z,t)
eta1 = Function('eta1')(x,y,z,t)
eta2 = Function('eta2')(x,y,z,t)

phi0 = Function('phi0')(x,y,z,t)
phi1 = Function('phi1')(x,y,z,t)
phi2 = Function('phi2')(x,y,z,t)

Mu0 = Function('Mu0')(x,y,z,t)
Mu1 = Function('Mu1')(x,y,z,t)
Mu2 = Function('Mu2')(x,y,z,t)

F0 = Function('F0')(x,y,z,t)
F1 = Function('F1')(x,y,z,t)
F2 = Function('F2')(x,y,z,t)


eta = eta0 + e*eta1 + e**2*eta2
phi = phi0 + e*phi1 + e**2*phi2
Mu  = Mu0 + e*Mu1 + e**2*Mu2
F   = F0 + e*F1 + e**2*F2

bodyEq = diff(Mu,t) - F
elevEq = diff(eta,t) - U*diff(eta,x) + diff(phi,x)*diff(eta,x) + diff(phi,y)*diff(eta,y) - diff(phi,z)
potEq  = diff(phi,t) - U*diff(phi,x) + 1/2*diff(phi,x)**2 + 1/2*diff(phi,y)**2 + 1/2*diff(phi,z)**2 + g*eta
lapEq  = diff(phi,x,2) + diff(phi,y,2) + diff(phi,z,2)

bodyEq0,bodyEq1 = expandEq(bodyEq)
elevEq0,elevEq1 = expandEq(elevEq)
potEq0,potEq1   = expandEq(potEq)
lapEq0,lapEq1   = expandEq(lapEq)

bodyEq0   = string_replacements(bodyEq0)
bodyEq1   = string_replacements(bodyEq1)
elevEq0   = string_replacements(elevEq0)
elevEq1   = string_replacements(elevEq1)
potEq0    = string_replacements(potEq0)
potEq1    = string_replacements(potEq1)
lapEq0    = string_replacements(lapEq0)
lapEq1    = string_replacements(lapEq1)


print(bodyEq0)
print(elevEq0)
print(potEq0)
print(lapEq0)
print("")
print(bodyEq1)
print(elevEq1)
print(potEq1)
print(lapEq1)
#
#
#print("r = ")
#print(rr)
#print("rx = ")
#print(rrx)
#print("ry = ")
#print(rry)
#print("rz = ")
#print(rrz)
#
#print("rxx = ")
#print(rrxx)
#print("ryy = ")
#print(rryy)
#print("rzz = ")
#print(rrzz)
#
#print("rxy = ")
#print(rrxy)
#print("rxz = ")
#print(rrxz)
#print("ryz = ")
#print(rryz)
#
#
#print("G = ")
#print(G)
#print("Gx = ")
#print(Gx)
#print("Gy = ")
#print(Gy)
#print("Gz = ")
#print(Gz)
#print("Gnx = ")
#print(Gnx)
#print("Gny = ")
#print(Gny)
#print("Gnz = ")
#print(Gnz)
#