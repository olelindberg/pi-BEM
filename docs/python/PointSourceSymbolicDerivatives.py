from sympy import *
from sympy.integrals.trigonometry import trigintegrate
from sympy.plotting import plot3d
from sympy.integrals.quadrature import gauss_legendre
import math

def string_replacements(strng):
    strng = str(strng).replace("r(x, y, z)","r")
    strng = str(strng).replace("Derivative","D")
    strng = str(strng).replace("D(r, x)","rx")
    strng = str(strng).replace("D(r, y)","ry")
    strng = str(strng).replace("D(r, z)","rz")
    strng = str(strng).replace("D(r, (x, 2))","rxx")
    strng = str(strng).replace("D(r, (y, 2))","ryy")
    strng = str(strng).replace("D(r, (z, 2))","rzz")
    strng = str(strng).replace("D(r, x, y)","rxy")
    strng = str(strng).replace("D(r, x, z)","rxz")
    strng = str(strng).replace("D(r, y, z)","ryz")
    strng = str(strng).replace("((x - x0)**2 + (y - y0)**2 + (z - z0)**2)**(-0.5)","1.0/r")
    strng = str(strng).replace("((x - x0)**2 + (y - y0)**2 + (z - z0)**2)**(-1.5)","1.0/r3")
    strng = str(strng).replace("1.0*","")
    strng = str(strng).replace("*1.0","")
    strng = str(strng).replace("(x - x0)","dx")
    strng = str(strng).replace("(y - y0)","dy")
    strng = str(strng).replace("(z - z0)","dz")
    strng = str(strng).replace("**2","2")
    strng = str(strng).replace("**3","3")
    strng = str(strng).replace("2*pi","two_pi")
    strng = str(strng).replace("4*pi","four_pi")
    return strng


nx,ny,nz,x,y,z,x0,y0,z0,pi = symbols('nx ny nz x y z x0 y0 z0 pi')
r = Function('r')(x,y,z)

rr    = ((x-x0)**2 + (y-y0)**2 + (z-z0)**2)**(1/2)

rrx    = diff(rr,x)
rry    = diff(rr,y)
rrz    = diff(rr,z)

rrxx    = diff(rr,x,x)
rryy    = diff(rr,y,y)
rrzz    = diff(rr,z,z)
rrxy    = diff(rr,x,y)
rrxz    = diff(rr,x,z)
rryz    = diff(rr,y,z)

G  = 1/(4*pi*r)

Gx    = diff(G,x)
Gy    = diff(G,y)
Gz    = diff(G,z)

Gn = nx*Gx + ny*Gy + nz*Gz

Gnx    = diff(Gn,x)
Gny    = diff(Gn,y)
Gnz    = diff(Gn,z)

rrx   = string_replacements(rrx)
rry   = string_replacements(rry)
rrz   = string_replacements(rrz)
rrxx   = string_replacements(rrxx)
rryy   = string_replacements(rryy)
rrzz   = string_replacements(rrzz)
rrxy   = string_replacements(rrxy)
rrxz   = string_replacements(rrxz)
rryz   = string_replacements(rryz)


G    = string_replacements(G)
Gx   = string_replacements(Gx)
Gy   = string_replacements(Gy)
Gz   = string_replacements(Gz)

Gnx   = string_replacements(Gnx)
Gny   = string_replacements(Gny)
Gnz   = string_replacements(Gnz)


print("r = ")
print(rr)
print("rx = ")
print(rrx)
print("ry = ")
print(rry)
print("rz = ")
print(rrz)

print("rxx = ")
print(rrxx)
print("ryy = ")
print(rryy)
print("rzz = ")
print(rrzz)

print("rxy = ")
print(rrxy)
print("rxz = ")
print(rrxz)
print("ryz = ")
print(rryz)


print("G = ")
print(G)
print("Gx = ")
print(Gx)
print("Gy = ")
print(Gy)
print("Gz = ")
print(Gz)
print("Gnx = ")
print(Gnx)
print("Gny = ")
print(Gny)
print("Gnz = ")
print(Gnz)
