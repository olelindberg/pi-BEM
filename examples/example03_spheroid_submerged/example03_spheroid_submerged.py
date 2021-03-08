import pyvista as pv
from pyvista import examples
import numpy as np


filename = "/home/ole/dev/projects/pi-BEM/build_cb/result_vector_results.vtu"
mesh = pv.read(filename)

print(mesh)
print(mesh.array_names)

# Make a geometric object to use as the glyph
geom = pv.Arrow()  # This could be any dataset

mesh.set_active_vectors('phi_gradient')
# Perform the glyph
glyphs = mesh.glyph(orient='phi_gradient', geom=geom)

# plot using the plotting class
p = pv.Plotter()
p.add_mesh(glyphs)
p.show()