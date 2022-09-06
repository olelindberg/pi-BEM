import numpy as np
from sympy import Naturals0
from HBIE_IntegralUtil import *

def HBIE_integral(mesh,shapeFunc,eta,n):

    derivative = mesh.derivatives(eta[0],eta[1])
    x_xi1_eta    = derivative[0]
    x_xi2_eta    = derivative[1]
    x_xi1xi1_eta = derivative[2]
    x_xi2xi1_eta = derivative[3]
    x_xi2xi2_eta = derivative[4]

    jac_eta     = mesh.jacobian(eta[0],eta[1])
    derivative  = mesh.jacobian_derivatives(eta[0],eta[1])
    jac_xi1_eta = derivative[0]
    jac_xi2_eta = derivative[1] 

    Na = shapeFunc.eval(eta)
    Na_xi1,Na_xi2 = shapeFunc.derivatives(eta)

    gaussx,gaussw = np.polynomial.legendre.leggauss(n)
    xi_vertices = np.array([[-1,1,1,-1],[-1,-1,1,1]])
    y = mesh.position(eta[0],eta[1])

    I0 = 0
    Im1 = 0
    Im2 = 0

    Im11 = 0
    Im12 = 0
    Im13 = 0
    Im14 = 0
    for i in [3,1,0,2]:

        v1 = xi_vertices[:,i]
        v2 = xi_vertices[:,(i+1)%4]
        if (i==3):

            print(v1,v2)
        r1 = v1 - eta
        r2 = v2 - eta

        theta1  = np.arctan2(r1[1],r1[0])
        theta2  = np.arctan2(r2[1],r2[0])
        if (theta2<theta1):
            theta2 += 2*np.pi

        for theta_gx,theta_gw in zip(gaussx,gaussw):

            s     = 0.5*(theta_gx+1)
            theta = (1-s)*theta1 + s*theta2
            if (i==3):
                print(theta)
            rho_hat = equation_of_external_contour_polar_coords(theta1,theta,eta,v1,v2)

            Ai = x_xi1_eta*np.cos(theta) + x_xi2_eta*np.sin(theta)
            Bi = 1/2*x_xi1xi1_eta*np.cos(theta)**2 + x_xi2xi1_eta*np.cos(theta)*np.sin(theta) + 1/2*x_xi2xi2_eta*np.sin(theta)**2
            if (i==3):
                print('x_xi1xi1_eta ',x_xi1xi1_eta)
                print('Ai, Bi       ',Ai,Bi)
            Ai = Ai.flatten()
            Bi = Bi.flatten()
            A = np.linalg.norm(Ai)
            C = np.inner(Ai,Bi)

            Ji0 = jac_eta
            Ji1 = jac_xi1_eta*np.cos(theta) + jac_xi2_eta*np.sin(theta)        
            if (i==3):
                print("Ji0,Ji1 ",Ji0,Ji1)
            Na0 = Na
            Na1 = Na_xi1*np.cos(theta) + Na_xi2*np.sin(theta)        

            beta = 1/A
            gamma = -C/A**4

            gi1 = Ai/A**2*(Bi@Ji0 + Ai@Ji1)

            bi0 = -Ji0
            bi1 = 3*gi1 - Ji1

            ai0 = np.outer(bi0,Na0)
            ai1 = np.outer(bi1,Na0) + np.outer(bi0,Na1)

            Sm3 = 1/A**3
            Sm2 = - 3*C/A**5

            Fm2 = -1/(4*np.pi)*Sm3*ai0
            Fm1 = -1/(4*np.pi)*(Sm2*ai0+Sm3*ai1)
            
            dtheta = (theta2-theta1)/2*theta_gw
            Im2 += - Fm2*(gamma/beta**2+1/rho_hat)*dtheta
            Im1 += Fm1*np.log(rho_hat/beta)*dtheta

            Im11 += Fm1*np.log(rho_hat/beta)
            Im12 += np.log(rho_hat/beta)
            Im13 += dtheta

            for rho_gx,rho_gw in zip(gaussx,gaussw):

                s   = 0.5*(rho_gx+1)
                rho = s*rho_hat

                xi1 = eta[0] + rho*np.cos(theta)
                xi2 = eta[1] + rho*np.sin(theta)
                xxx = mesh.position(xi1,xi2)

                rvec = xxx-y
                r = np.linalg.norm(rvec)
                ri = rvec/r
                Ji = mesh.jacobian(xi1,xi2)
                J  = np.linalg.norm(Ji)
                ni = Ji/J
                N = shapeFunc.eval(np.array([xi1,xi2]))
                F = - 1/(4*np.pi*r**3)*np.outer(3*ri*np.inner(ri,ni) - ni,N)*J*rho

                drho = rho_hat/2*rho_gw
                I0 += (F - (Fm2/rho**2 + Fm1/rho))*drho*dtheta
#        print(Im1[2])
        #print(" ")
        if (i==3):
            print("Im11 ", Im11[2][0])
        #print(Im12)

#    print(I0)
#    print(Im2)

    I = I0+Im1+Im2
    return I
