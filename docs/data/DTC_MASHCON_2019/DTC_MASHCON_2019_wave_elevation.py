import numpy as np
import matplotlib.pyplot as plt

savefigures = False

elev = np.genfromtxt("C3/output/elevation.csv",delimiter=",")

g     = 9.80665
Fn    = 0.25
L     = 2.5
h     = 2.0
U     = Fn*np.sqrt(g*L)
Fnh     = U/np.sqrt(g*h)

plt.plot(elev[:,0], -elev[:,2],"b.")
#plt.title(r'Wave elevation, Wigley Hull, ITTC 1984, $Fn_L$=' + "{:1.2f}".format(Fn) + r', $Fn_h$=' + "{:1.2f}".format(Fnh))
#plt.xlabel(r'$x/L$')
#plt.ylabel(r'$\eta/(u^2/g)$')
plt.grid(True)

#if savefigures:
#    filename = 'wigleyhull_wave_elevation_ITTC_1984.pdf'
#    plt.savefig(filename,bbox_inches='tight')
#    filename = 'wigleyhull_wave_elevation_ITTC_1984.png'
#    plt.savefig(filename,bbox_inches='tight')

plt.show()