import pyvista as pv
import numpy as np

show = True
showPotential = True
show_edges = False
zoom = 0.29

edge_color = 'black'
density = 1000
velocity = 4.499466057211677
filenames = ["/home/ole/dev/projects/pi-BEM/docs/data/KVLCC2_SHAL_FHR_SIMMAN_2014/T0Z02A03_C4/mesh07/output/result_scalar_results.vtu"]

sargs = dict(
    title_font_size=20,
    label_font_size=30,
    color="black",
    shadow=False,
    n_labels=5,
    italic=False,
    fmt="%.2f",
    font_family="arial",
    position_x=0.05, position_y=0.05
)

stagnationPressure = 0.5*density*velocity**2


for filename in filenames:

    print(filename)
    mesh = pv.read(filename)

    if show:
        p = pv.Plotter()
    else:
        p = pv.Plotter(off_screen=True)

    p.set_background("white")
    p.camera_position = [(zoom*1000, zoom * 200, zoom * 150),
                         (180, 0, 0), (0, 0, 1)]

    if showPotential:
        press = mesh.get_array('hydrodynamic_pressure')/stagnationPressure
        p.add_mesh(mesh, scalars=press,
                   show_edges=show_edges, edge_color=edge_color, smooth_shading=True, opacity=1, cmap='jet', clim=(-1, 1), scalar_bar_args=sargs)
    else:
        p.add_mesh(mesh, color='white',
                   show_edges=show_edges, edge_color='black', opacity=1)

    if show:
        p.show()
    else:

        imagename = 'KCS_SHAL_FHR_SIMMAN_2014_Edges' + str(show_edges) + \
            'Potential' + str(showPotential) + ".png"
        print(imagename)
        p.screenshot(imagename, window_size=(1920, 1080))
        p.clear()
