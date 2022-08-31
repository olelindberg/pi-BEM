from sympy import *

rho, A, B, C, Ai, Aj, Bi, Bj, Jj0, Jj1, Ji0, Ji1,N0c,N1c = symbols("rho A B C Ai Aj Bi Bj Jj0 Jj1 Ji0 Ji1 N0c N1c")

#1/(4*pi)*
Fi = -(1/(rho**3*A**3) - 3*C/(rho**2*A**5) + 1/rho)*((Ai/A + rho*(Bi/A-Ai*C/A**3) + rho**2)*(Aj/A + rho*Bj/A-Aj*C/A**3 + rho**2)*(Jj0 + rho*Jj1 + rho**2) -(Ji0 + rho*Ji1 + rho**2))*(N0c + rho*N1c + rho**2)*rho

Fi = expand(Fi)

for i in range(-2,0):
    print(i)
    print(collect(Fi.coeff(rho,i),A))
