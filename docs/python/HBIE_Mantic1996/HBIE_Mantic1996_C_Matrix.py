from sympy import *
init_printing(use_unicode=True)
r, phi, theta = symbols("r phi theta")

x = r*cos(phi)*sin(theta)
y = r*sin(phi)*sin(theta)
z = r*cos(theta)
       
jac = Matrix([  [diff(x,r),diff(x,theta),diff(x,phi)],
                [diff(y,r),diff(y,theta),diff(y,phi)],
                [diff(z,r),diff(z,theta),diff(z,phi)]]).det()

xx = Matrix([[x*x,x*y,x*z],
            [y*x,y*y,y*z],
            [z*x,z*y,z*z]])

C11 = 3/(4*pi)*integrate(integrate(integrate(xx*jac,(phi,0,3/2*pi)),(theta,0,pi)),(r,0,1))



print(C11)