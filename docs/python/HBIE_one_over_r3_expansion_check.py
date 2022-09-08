from cProfile import label
import numpy as np
import matplotlib.pyplot as plt


data = np.array([[4.2612e-03,   8.4159e-03,   1.6776e+06,   1.6776e+06,  -8.5441e+00],
                 [2.2342e-02,   4.3821e-02,   1.1884e+04,   1.1886e+04,  -1.6919e+00],
                 [5.4427e-02,   1.0550e-01,   8.5173e+02,   8.5247e+02,  -7.3814e-01],
                 [9.9770e-02,   1.9039e-01,   1.4491e+02,   1.4534e+02,  -4.3412e-01],
                 [1.5731e-01,   2.9480e-01,   3.9030e+01,   3.9328e+01,  -2.9818e-01],
   [2.2570e-01,   4.1494e-01,   1.3997e+01,   1.4222e+01,  -2.2420e-01],
   [3.0333e-01,   5.4719e-01,   6.1035e+00,   6.2816e+00,  -1.7814e-01],
   [3.8839e-01,   6.8829e-01,   3.0668e+00,   3.2133e+00,  -1.4655e-01],
   [4.7888e-01,   8.3522e-01,   1.7163e+00,   1.8396e+00,  -1.2335e-01],
   [5.7268e-01,   9.8509e-01,   1.0461e+00,   1.1517e+00,  -1.0556e-01],
   [6.6760e-01,   1.1350e+00,   6.8392e-01,   7.7551e-01,  -9.1591e-02],
   [7.6140e-01,   1.2820e+00,   4.7456e-01,   5.5509e-01,  -8.0527e-02],
   [8.5189e-01,   1.4233e+00,   3.4682e-01,   4.1860e-01,  -7.1777e-02],
   [9.3695e-01,   1.5560e+00,   2.6544e-01,   3.3036e-01,  -6.4921e-02],
   [1.0146e+00,   1.6774e+00,   2.1186e-01,   2.7149e-01,  -5.9630e-02],
   [1.0830e+00,   1.7851e+00,   1.7580e-01,   2.3143e-01,  -5.5633e-02],
   [1.1405e+00,   1.8765e+00,   1.5134e-01,   2.0405e-01,  -5.2703e-02],
   [1.1859e+00,   1.9493e+00,   1.3501e-01,   1.8566e-01,  -5.0651e-02],
   [1.2179e+00,   2.0014e+00,   1.2474e-01,   1.7407e-01,  -4.9328e-02],
   [1.2360e+00,   2.0309e+00,   1.1937e-01,   1.6800e-01,  -4.8627e-02]])

rho                 = data[:,0]
r                   = data[:,1]
r3_inv              = data[:,2]
r3_inv_expansion    = data[:,3]
err                 = data[:,4]


plt.plot(np.log10(rho),np.log10(r3_inv),'r',label=r"$\log_{10}(1/r^3)$")
plt.plot(np.log10(rho),np.log10(r3_inv_expansion),'b',label=r"$\log_{10}(1/r^3), expansion$")
plt.plot(np.log10(rho),np.log10(np.fabs(err)),'g',label=r"$\log_{10}(|err|)$")
plt.plot(np.log10(rho),np.log10(1/rho),'k--',label=r"$\log_{10}(1/\rho)$")
plt.legend()
plt.grid(True)
plt.xlabel(r"$\log_{10}(\rho)$")

plt.show()