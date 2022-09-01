from sympy import *

xi1,xi2 = symbols("xi1 xi2")

phi     = pi/4*(xi1+1)/2
theta   = pi/2-pi/4*(xi2+1)/2

x = cos(phi)*sin(theta) 
y = sin(phi)*sin(theta)
z = cos(phi)*cos(theta) 

x_xi1    = [x.diff(xi1),y.diff(xi1),z.diff(xi1)]
x_xi2    = [x.diff(xi2),y.diff(xi2),z.diff(xi2)]
x_xi1xi1 = [x.diff(xi1,2),y.diff(xi1,2),z.diff(xi1,2)]
x_xi1xi2 = [x.diff(xi1).diff(xi2),y.diff(xi1).diff(xi2),z.diff(xi1).diff(xi2)]
x_xi2xi2 = [x.diff(xi2,2),y.diff(xi2,2),z.diff(xi2,2)]

jac = [x_xi1[1]*x_xi2[2] - x_xi1[2]*x_xi2[1],       
       x_xi1[2]*x_xi2[0] - x_xi1[0]*x_xi2[2],
       x_xi1[0]*x_xi2[1] - x_xi1[1]*x_xi2[0]]
jac_xi1 = [jac[0].diff(xi1),jac[1].diff(xi1),jac[2].diff(xi1)]
jac_xi2 = [jac[0].diff(xi2),jac[1].diff(xi2),jac[2].diff(xi2)]

print("x_xi1    = ",x_xi1   ) 
print("x_xi2    = ",x_xi2   ) 
print("x_xi1xi1 = ",x_xi1xi1) 
print("x_xi1xi2 = ",x_xi1xi2) 
print("x_xi2xi2 = ",x_xi2xi2) 

print("jac     = ",jac) 
print("jac_xi1 = ",jac_xi1) 
print("jac_xi2 = ",jac_xi2) 