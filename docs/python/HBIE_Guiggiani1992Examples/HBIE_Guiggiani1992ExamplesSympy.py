from sympy import *

xi1,xi2 = symbols("xi1 xi2")

angle = pi/2*(xi2+1)/2
x = cos(angle) 
y = xi1 + 1
z = sin(angle) 

x_xi1    = x.diff(xi1)
y_xi1    = y.diff(xi1)
z_xi1    = z.diff(xi1)
   
x_xi2    = x.diff(xi2)
y_xi2    = y.diff(xi2)
z_xi2    = z.diff(xi2)

x_xi1xi1 = x.diff(xi1,2)
y_xi1xi1 = y.diff(xi1,2)
z_xi1xi1 = z.diff(xi1,2)

x_xi2xi2 = x.diff(xi2,2)
y_xi2xi2 = y.diff(xi2,2)
z_xi2xi2 = z.diff(xi2,2)

x_xi1xi2 = x.diff(xi1).diff(xi2)
y_xi1xi2 = y.diff(xi1).diff(xi2)
z_xi1xi2 = z.diff(xi1).diff(xi2)


jac = [y_xi1*z_xi2 - z_xi1*y_xi2,       z_xi1*x_xi2 - x_xi1*z_xi2,       x_xi1*y_xi2 - y_xi1*x_xi2]
jac_xi1 = [jac[0].diff(xi1),jac[1].diff(xi1),jac[2].diff(xi1)]
jac_xi2 = [jac[0].diff(xi2),jac[1].diff(xi2),jac[2].diff(xi2)]



print("x_xi1    = ",x_xi1   ) 
print("y_xi1    = ",y_xi1   ) 
print("z_xi1    = ",z_xi1   ) 
print("x_xi2    = ",x_xi2   ) 
print("y_xi2    = ",y_xi2   ) 
print("z_xi2    = ",z_xi2   ) 
print("x_xi1xi1 = ",x_xi1xi1) 
print("y_xi1xi1 = ",y_xi1xi1) 
print("z_xi1xi1 = ",z_xi1xi1) 
print("x_xi2xi2 = ",x_xi2xi2) 
print("y_xi2xi2 = ",y_xi2xi2) 
print("z_xi2xi2 = ",z_xi2xi2) 
print("x_xi1xi2 = ",x_xi1xi2) 
print("y_xi1xi2 = ",y_xi1xi2) 
print("z_xi1xi2 = ",z_xi1xi2) 

print("jac     = ",jac) 
print("jac_xi1 = ",jac_xi1) 
print("jac_xi2 = ",jac_xi2) 