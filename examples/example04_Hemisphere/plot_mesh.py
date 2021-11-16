import pyvista as pv


filename = "output/result_scalar_results.vtu"
mesh = pv.read(filename)
#cpos = mesh.plot(show_edges=True)


centers = mesh.cell_centers()

p = pv.Plotter()
p.add_mesh(mesh, show_edges=True, line_width=1)
#p.add_mesh(centers, color="r", point_size=4.0, render_points_as_spheres=True)

point_labels = [f"{i}" for i in range(centers.n_points)]

p.add_point_labels(centers, labels=point_labels, shape_opacity=0,
                   bold=False, render_points_as_spheres=True, point_color='r')
p.show()
