import numpy as np
import matplotlib.pyplot as plt

savefigures = True

elev = np.genfromtxt("example14_wigleyhull_ITTC_1984/output/elevation.csv",delimiter=",")

g     = 9.80665
Fn    = 0.25
L     = 2.5
U     = Fn*np.sqrt(g*L)
x_IHI = np.array([-1,-0.95,-0.9,-0.85,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.95,1])
h_IHI = np.array([0.177,0.378,0.374,0.281,0.151,-0.077,-0.174,-0.191,-0.120,-0.050,-0.016,0.022,0.014,-0.027,-0.073,-0.086,-0.080,-0.056,-0.046,-0.003,0.012,0.013,0.027,0.070,0.187])


plt.plot((x_IHI+1)/2, h_IHI,"k-")
plt.plot(elev[:,0]/2.5, g/(U*U)*elev[:,2],"b.")
plt.xlabel(r'$x/L$')
plt.ylabel(r'$\eta/(u^2/g)$')
plt.grid(True)

if savefigures:
    filename = 'example14_wigleyhull_ITTC_1984/wigleyhull_wave_elevation_ITTC_1984.pdf'
    plt.savefig(filename,bbox_inches='tight')
    filename = 'example14_wigleyhull_ITTC_1984/wigleyhull_wave_elevation_ITTC_1984.png'
    plt.savefig(filename,bbox_inches='tight')

plt.show()