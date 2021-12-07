import pyvista as pv
import numpy as np

showCellCenters = True

filepath = "/home/ole/dev/projects/pi-BEM/build"

rho = 1000
V = 3.0830504374726018
pressureScale = 0.5*rho*V*V

numProc = 4
pid = 3
filename = "/home/ole/dev/projects/pi-BEM/docs/data/BankEffects_JMST_2019/case1/mesh5/output/result_scalar_results.vtu"

mesh = pv.read(filename)
print(mesh.array_names)

pres = mesh.get_array('hydrodynamic_pressure')/pressureScale

p = pv.Plotter()
p.add_mesh(mesh, scalars=pres, cmap='jet',
           clim=[-1, 1], smooth_shading=True, n_colors=32, ambient=0.3)

zoom = 0.6
p.camera_position = [(zoom*1000, zoom * 300, zoom * 320),
                     (200, 0, 0), (0, 0, 1)]

p.show(screenshot='BankEffectsForKVLCC2_case1.png')
