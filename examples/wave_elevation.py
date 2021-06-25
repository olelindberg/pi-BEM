import numpy as np
import matplotlib.pyplot as plt

savefigures = False

elev = np.genfromtxt("example05_DTC_ship/output/elevation.csv",delimiter=",")
#elev = np.genfromtxt("example14_wigleyhull_ITTC_1984/output/elevation.csv",delimiter=",")

g     = 9.80665
Fn    = 0.1
L     = 355.0
U     = Fn*np.sqrt(g*L)

#plt.plot((x_IHI+1)/2, h_IHI,"k-")
plt.plot(elev[:,0]/L, g/(U*U)*elev[:,2],"b.")
plt.xlabel(r'$x/L$')
plt.ylabel(r'$\eta/(u^2/g)$')
plt.grid(True)

if savefigures:
    filename = 'example14_wigleyhull_ITTC_1984/wigleyhull_wave_elevation_ITTC_1984.pdf'
    plt.savefig(filename,bbox_inches='tight')
    filename = 'example14_wigleyhull_ITTC_1984/wigleyhull_wave_elevation_ITTC_1984.png'
    plt.savefig(filename,bbox_inches='tight')

plt.show()