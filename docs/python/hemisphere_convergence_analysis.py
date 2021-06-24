import pyvista as pv
import math
import numpy as np
import matplotlib.pyplot as plt

r  = 1
U  = -1
mu = 2*math.pi*U*r**3

hgrid = []
err_max = []
err_rms = []

for i in range(1,7):

  print(i)

  filename  = '/home/ole/dev/projects/pi-BEM/docs/data/hemisphere0' + str(i) + '/output/result_scalar_results.vtu'
  mesh      = pv.read(filename)

  mesh = mesh.compute_cell_quality(quality_measure='area')
  hgrid.append(np.sqrt(np.max(mesh.cell_arrays['CellQuality'])))

  x         = mesh.points[:,0]
  y         = mesh.points[:,1]
  z         = mesh.points[:,2]

  phi       = mesh.point_arrays['phi']
  phi_ex    = mu*x/(4*math.pi*(x**2 + y**2 + z**2)**(3/2))

  offset = phi_ex[0]-phi[0]

  phi = phi+offset

  err = phi_ex - phi
  err_rms.append(np.sqrt(np.mean(err**2)))
  err_max.append(np.max(np.abs(err)))

fit_rms = np.polyfit(np.log10(hgrid),np.log10(err_rms), 1)
fit_max = np.polyfit(np.log10(hgrid),np.log10(err_max), 1)
print(fit_max)
print(fit_rms)
plt.plot(np.log10(hgrid),np.log10(err_rms),'r-o',label='RMS')
plt.plot(np.log10(hgrid),np.log10(err_max),'b-o',label='INF')
plt.plot(np.log10(hgrid),np.log10(hgrid)-0.5,'k--',label='slope = 1')
plt.plot(np.log10(hgrid),2*np.log10(hgrid)-0.75,'k-.',label='slope = 2')
plt.grid(True)
plt.legend()
plt.xlabel(r'$log_{10}(h)$')
plt.ylabel(r'$log_{10}(err)$')
plt.savefig('hemisphere_convergence_analysis.pdf')
plt.show()
  # Plot stuff:
#  plotter = pv.Plotter()    # instantiate the plotter
#  plotter.add_mesh(mesh,scalars=err,show_edges=True,colormap="jet")    # add a mesh to the scene
#  plotter.add_scalar_bar()
#  cpos = plotter.show()

