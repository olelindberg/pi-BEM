import pyvista as pv
import numpy as np


filepath = "/home/ole/dev/projects/pi-BEM/build"

numProc = 2
pid = 1
filename = "scalars_np_" + str(numProc) + "_pid_" + str(pid) + ".vtu"

print(filepath)
print(filename)
mesh = pv.read(filepath + "/" + filename)

filename2 = "cellcenters_np_" + str(numProc) + "_pid_" + str(pid) + ".csv"
data = np.genfromtxt(filepath + "/" + filename2, delimiter=' ')
print(data)

centers = mesh.cell_centers()
point_labels = [f"{i}" for i in range(centers.n_points)]

cellcenters = pv.PolyData()
cellcenters.points = data[:, 1:4]
point_labels = [f"{int(i)}" for i in data[:, 0]]

centers = cellcenters

p = pv.Plotter()
p.add_mesh(mesh, show_edges=True, line_width=1, color='w')
p.add_mesh(centers, color="r", point_size=4.0,
           render_points_as_spheres=True)


p.add_point_labels(centers, labels=point_labels, shape_opacity=0,
                   bold=False, render_points_as_spheres=True, point_color='r')
p.add_title(filename)
p.show()
