from sympy import *
r, phi, theta = symbols("r phi theta")

phi_min = 0
phi_max = 3/2*pi

theta_min = 0
theta_max = pi

# Spherical coordinates:
x = r*cos(phi)*sin(theta)
y = r*sin(phi)*sin(theta)
z = r*cos(theta)

# Jacobian of the coordinate transform:       
jac = Matrix([  [diff(x,r),diff(x,theta),diff(x,phi)],
                [diff(y,r),diff(y,theta),diff(y,phi)],
                [diff(z,r),diff(z,theta),diff(z,phi)]]).det()

# Outer product:
xx = Matrix([[x*x,x*y,x*z],
             [y*x,y*y,y*z],
             [z*x,z*y,z*z]])

# Integrate over unit sphere:
C = 3/(4*pi)*integrate(integrate(xx*jac,(phi,phi_min,phi_max)),(theta,theta_min,theta_max))

# Insert unit radius r=1:
C = C.subs(r,1)

# Show result:
print(C)