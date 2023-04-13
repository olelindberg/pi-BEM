from sympy import *

# Derivation of the double body equations consist of the following steps:
# 1) Moving frame of ref.
# 2) Linearization around still water level (SWL), z=z_SWL.
# 3) Pertubation expansion.


def string_replacements(strng):
    strng = str(strng).replace("(x, y, z, t)","")
    strng = str(strng).replace("Derivative","D")
    strng = str(strng).replace("1.0*","")

    return strng

def ExpandEquation(eq):

    eq = expand(eq)
    eq = collect(eq,e)
    eq0 = eq.coeff(e,0)
    eq1 = eq.coeff(e,1)

    return eq0,eq1

nx,ny,nz,x,y,z,x0,y0,z0,pi,t,e,U,g = symbols('nx ny nz x y z x0 y0 z0 pi t e U g')

eta = Function('eta')(x,y,t)
eta0 = Function('eta0')(x,y,t)
eta1 = Function('eta1')(x,y,t)

phi = Function('phi')(x,y,z,t)
phi0 = Function('phi0')(x,y,z,t)
phi1 = Function('phi1')(x,y,z,t)



phi_taylor_exp = phi# + eta*diff(phi,z)

print("Taylor exp of phi around z_SWL:")
print(phi_taylor_exp)

bernoulli_eq  = diff(phi,t) - U*diff(phi,x) + 1/2*diff(phi,x)**2 + 1/2*diff(phi,y)**2 + 1/2*diff(phi,z)**2 + g*eta

print("Bernoulli eq:")
print(bernoulli_eq)

bernoulli_eq = bernoulli_eq.subs(phi,phi_taylor_exp)

print("Bernoulli after subs with Taylor exp:")
print(bernoulli_eq)

eta_pertubation_exp = eta0# + e*eta1
phi_pertubation_exp = phi0# + e*phi1

bernoulli_eq = bernoulli_eq.subs(phi,phi_pertubation_exp)
bernoulli_eq = bernoulli_eq.subs(eta,eta_pertubation_exp)

print("Bernoulli after subs with pertubation exp:")
print(bernoulli_eq)


print("Bernoulli O(e^0) terms:")
print(bernoulli_eq.expand())
#print("Bernoulli O(e^1) terms:")
#print(bernoulli_eq1)



# elevEq = diff(eta,t) - U*diff(eta,x) + diff(phi,x)*diff(eta,x) + diff(phi,y)*diff(eta,y) - diff(phi,z)
# lapEq  = diff(phi,x,2) + diff(phi,y,2) + diff(phi,z,2)

# elevEq0,elevEq1 = expandEq(elevEq)
# lapEq0,lapEq1   = expandEq(lapEq)

# elevEq0   = string_replacements(elevEq0)
# elevEq1   = string_replacements(elevEq1)
# lapEq0    = string_replacements(lapEq0)
# lapEq1    = string_replacements(lapEq1)


# print(elevEq0)
# print(lapEq0)
# print("")
# print(elevEq1)
# print(lapEq1)
