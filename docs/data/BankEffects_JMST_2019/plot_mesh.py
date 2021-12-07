import pyvista as pv
import numpy as np

showCellCenters = True

filepath = "/home/ole/dev/projects/pi-BEM/build"

numProc = 4
pid = 3
filename = "/home/ole/dev/projects/pi-BEM/docs/data/BankEffects_JMST_2019/case3/mesh5/output/result_scalar_results.vtu"

mesh = pv.read(filename)


centers = mesh.cell_centers()
point_labels = [f"{i}" for i in range(centers.n_points)]

p = pv.Plotter()

p.add_mesh(mesh, show_edges=True, line_width=1)
p.add_mesh(centers, color="r", point_size=4.0,
           render_points_as_spheres=True)

p.add_point_labels(centers, labels=point_labels, shape_opacity=0,
                   bold=False, render_points_as_spheres=True, point_color='r')

# p.add_point_labels(nodes, labels=nodes_labels, shape_opacity=0,
#                   bold=False, render_points_as_spheres=True, point_color='b')

p.add_title(filename)
p.show()
