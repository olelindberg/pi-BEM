import pyvista as pv
import numpy as np
import glob
import os

from pyvista.plotting import scalar_bars

numsteps = 1

show = False
showPotential = True
show_edges = True

zoomType = 0
zoom1 = 0.3
zoom2 = 1.0

edge_color = 'black'

density = 1000
velocity = 4.499466057211677

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

filenames = "/home/ole/dev/projects/pi-BEM/docs/data/KCS_SHAL_FHR_SIMMAN_2014/C0101A03_C6_movie/output/result_scalar_results.vtu"
#filenames = "/home/ole/dev/temp/scalars*.vtu"
#filenames = "/home/ole/dev/temp/trimesh*.vtu"
#filenames = "/home/ole/dev/temp/trimesh000002.vtu"

dirs = sorted(glob.glob(filenames))
print(dirs)
dirs2 = []
for ii in range(numsteps):
    dirs2 = dirs2 + dirs
print(dirs2)


cnt = 0.0
for filename in dirs2:
    print(filename)
    mesh = pv.read(filename)

    if show:
        p = pv.Plotter()
    else:
        p = pv.Plotter(off_screen=True)

    p.set_background("white")
    if showPotential:
        press = mesh.get_array('hydrodynamic_pressure')/stagnationPressure

        p.add_mesh(mesh, scalars=press,
                   show_edges=show_edges, edge_color=edge_color, smooth_shading=True, opacity=1, cmap='jet', clim=(-1, 1), scalar_bar_args=sargs)
        print(p.scalar_bars.values())
    else:
        p.add_mesh(mesh, color='white',
                   show_edges=show_edges, edge_color='black', opacity=1)
    wgt = 0
    if len(dirs2) > 1:
        wgt = cnt/(len(dirs2)-1)
    zoom = zoom1*(1.0-wgt) + zoom2*wgt

    p.camera_position = [(zoom*1000, zoom * 400, zoom * 300),
                         (155, 0, 0), (0, 0, 1)]
    if show:
        p.show()
    else:

        cntstr = ''
        if numsteps > 1:
            cntstr = '{0:06d}'.format(int(cnt))

        imagename = os.path.splitext(filename)[0] + cntstr + ".png"
        print(imagename)
        p.screenshot(imagename, window_size=(1920, 1080))
        p.clear()

    cnt = cnt+1
