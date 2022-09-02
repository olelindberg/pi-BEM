import numpy as np

def LagrangePolynomial2DDerivativeMatrices(x,derivative):

    D = LagrangePolynomialMthDerivativeMatrix(x,derivative)
    [D1,D2] = LagrangePolynomial2DDerivativeMatricesAssemble(D)
    return D1,D2

def LagrangePolynomial2DDerivativeMatricesAssemble(D):

    N1D = D.shape[0]
    N2D = N1D*N1D

    # Make Dxi1 matrix:
    Dxi1 = np.zeros((N2D,N2D))
    for k in range(N1D):
        id = np.arange(k*N1D,(k+1)*N1D,1)
        for j in range(N1D):
            for i in range(N1D):
                Dxi1[id[i],id[j]] = D[i,j]

    # Make Dxi2 matrix:
    Dxi2 = np.zeros((N2D,N2D))
    for k in range(N1D):
        for j in range(N1D):
            for i in range(N1D): # Nodes loop:
                Dxi2[j*N1D+k,i*N1D+k] = D[j,i]

    return Dxi1,Dxi2

def LagrangePolynomial1DBarycentricWeights(xi):

    #--------------------------------------------------------------------------
    # Barycentric weights (see Kopriva 2009 alg. 30):
    #--------------------------------------------------------------------------
    bw    = 1+0*xi
    order = len(xi)-1

    #----------------------------------------------------------------------------
    # Calculate the weights:
    #----------------------------------------------------------------------------
    for i in range(1,order+1):
        for j in range(i):
            bw[j] = bw[j]*(xi[j] - xi[i])
            bw[i] = bw[i]*(xi[i] - xi[j])

    #----------------------------------------------------------------------------
    # Normalize the weights:
    #----------------------------------------------------------------------------
    bw = 1.0/bw

    return bw

def LagrangePolynomial1DInterpMatrix(xi,xi0):

    #------------------------------------------------------------------------------
    # Lagrange interpolating polynomials matrix (see Kopriva 2009 alg. 34):
    #------------------------------------------------------------------------------

    order = len(xi)-1

    bw = LagrangePolynomial1DBarycentricWeights(xi)
    LP = np.zeros((len(xi0),len(xi)))

    for j in range(len(xi0)):

        lp          = 0*xi
        node_match  = False
        for i in range(order+1):
            lp[i] = 0.0
            if (np.fabs(xi[i]-xi0[j])<1e-12):
                lp[i]      = 1.0
                node_match = True

        if (not node_match):
            for i in range(order+1):
                lp[i]   = bw[i]/(xi0[j]-xi[i])
            
            lp = lp/sum(lp)
        

        LP[j,:] = lp
    return LP

def LagrangePolynomialDerivativeMatrix(xi):

    #------------------------------------------------------------------------------
    # Polynomial derivative matrix (see Kopriva 2009 alg. 37):
    #------------------------------------------------------------------------------

    # Assumes: 1) node positions, xi.
    #          2) bary_centric_weights, bw.

    bw = LagrangePolynomial1DBarycentricWeights(xi)

    # Allocate derivative matrix:
    order   = len(xi) - 1
    D       = np.zeros((order+1,order+1))

    # Calculate derivative matrix:
    for i in range(order+1):
        
        D[i,i] = 0.0
        
        for j in range(order+1):
            if i != j:
                
                D[i,j] = bw[j]/bw[i]*1.0/(xi[i]-xi[j])
                
                # Negative sum trick (see Kopriva 2009 page 55):
                D[i,i] = D[i,i] - D[i,j]
    return D,bw


def LagrangePolynomial2DInterpMatrix(xi1,xi2,xi10,xi20):

    interp_xi1 = LagrangePolynomial1DInterpMatrix(xi1,np.array([xi10]))
    interp_xi2 = LagrangePolynomial1DInterpMatrix(xi2,np.array([xi20]))
    
    N1 = len(xi1)
    N2 = len(xi2)

    # Make Interp1 matrix:
    interp1 = np.zeros((N2,N1*N2))
    for j in range(N1):
        for i in range(N2):
            interp1[i,j+i*N1] = interp_xi1.flatten()[j]

    interp2d = interp_xi2@interp1
    return interp2d

def LagrangePolynomialMthDerivativeMatrix(xi,M):

#------------------------------------------------------------------------------
# Polynomial derivative matrix (see Kopriva 2009 alg. 38):
#------------------------------------------------------------------------------
# Assumes:  1) node positions, xi.
#           2) bary_centric_weights, bw.
#           3) Lagrange polynomial derivative matrix

    D,bw = LagrangePolynomialDerivativeMatrix(xi)
    N  = len(xi) - 1

    for k in range(2,M+1):
        Dk = 0*D
        for i in range(N):
            for j in range(N):
                if i != j:
                    Dk[i,j] = k/(xi[i]-xi[j])*(bw[j]/bw[i]*D[i,i] - D[i,j])
                    # Negative sum trick (see Kopriva 2009 page 55):
                    Dk[i,i] = Dk[i,i] - Dk[i,j]
        D = Dk
    return D

