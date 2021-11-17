import pyvista as pv
import numpy as np

showCellCenters = True

filepath = "/home/ole/dev/projects/pi-BEM/build"

numProc = 2
pid = 1
filename = "scalars_np_" + str(numProc) + "_pid_" + str(pid) + ".vtu"

print(filepath)
print(filename)
mesh = pv.read(filepath + "/" + filename)

filename2 = "cellcenters_np_" + str(numProc) + "_pid_" + str(pid) + ".csv"
data2 = np.genfromtxt(filepath + "/" + filename2, delimiter=' ')

filename3 = "nodes_np_" + str(numProc) + "_pid_" + str(pid) + ".csv"
data3 = np.genfromtxt(filepath + "/" + filename3, delimiter=' ')


#centers = mesh.cell_centers()
#point_labels = [f"{i}" for i in range(centers.n_points)]

cellcenters = pv.PolyData()
cellcenters.points = data2[:, 1:4]
cellcenter_labels = [f"{int(i)}" for i in data2[:, 0]]

cellcenter = cellcenters

nodes = pv.PolyData()
nodes.points = data3[:, 1:4]
nodes_labels = [f"{int(i)}" for i in data3[:, 0]]



p = pv.Plotter()

p.add_mesh(mesh, show_edges=True, line_width=1)

p.add_mesh(cellcenter, color="r", point_size=4.0,
           render_points_as_spheres=True)

if showCellCenters:
    p.add_point_labels(cellcenter, labels=cellcenter_labels, shape_opacity=0,
                    bold=False, render_points_as_spheres=True, point_color='r')

p.add_point_labels(nodes, labels=nodes_labels, shape_opacity=0,
                   bold=False, render_points_as_spheres=True, point_color='b')

p.add_title(filename)
p.show()
