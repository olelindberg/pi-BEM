import pyvista as pv
from pyvista import examples
import numpy as np
import math
from matplotlib import cm

def FlowAroundSphere(x,y,z,R,Ubody,Uinf):
    U   = Uinf - Ubody
    mu  = 2*math.pi*U*R**3
    phi = Uinf*x + mu/(4*math.pi)*x*(x**2 + y**2 + z**2)**(-3/2)
    u   = Uinf - 0.75*mu*x**2*(x**2 + y**2 + z**2)**(-2.5)/math.pi + mu*(x**2 + y**2 + z**2)**(-1.5)/(4*math.pi)
    v   = -0.75*mu*x*y*(x**2 + y**2 + z**2)**(-2.5)/math.pi
    w   = -0.75*mu*x*z*(x**2 + y**2 + z**2)**(-2.5)/math.pi
    return [phi,u,v,w]

R     = 1
Uinf  = 0
Ubody = 1

colormap = cm.get_cmap('jet')

filename    = "/home/ole/dev/projects/pi-BEM/build_cb/result_scalar_results.vtu"
mesh_scalar = pv.read(filename)
print(mesh_scalar)
print(mesh_scalar.array_names)


filename    = "/home/ole/dev/projects/pi-BEM/build_cb/result_vector_results.vtu"
mesh_vector = pv.read(filename)
print(mesh_vector)
print(mesh_vector.array_names)

phi       = mesh_scalar.get_array(mesh_scalar.array_names[0])
vel       = mesh_vector.get_array(mesh_vector.array_names[0])
phi_exact = phi*0
vel_exact = vel*0

cnt = 0
for point in mesh_scalar.points:
    state = FlowAroundSphere(point[0],point[1],point[2],R,Ubody,Uinf)
    phi_exact[cnt]   = state[0]
    vel_exact[cnt,0] = state[1]
    vel_exact[cnt,1] = state[2]
    vel_exact[cnt,2] = state[3]
    cnt = cnt + 1

mesh_scalar.add_field_array(phi_exact,        name="phi_exact")
mesh_scalar.add_field_array(vel_exact[:,0],   name="u_exact")
mesh_scalar.add_field_array(vel_exact[:,1],   name="v_exact")
mesh_scalar.add_field_array(vel_exact[:,2],   name="w_exact")
mesh_scalar.add_field_array(vel[:,0],         name="u")
mesh_scalar.add_field_array(vel[:,1],         name="v")
mesh_scalar.add_field_array(vel[:,2],         name="w")
mesh_scalar.add_field_array(vel,              name="vel")
mesh_scalar.add_field_array(vel_exact,        name="vel_exact")

print(mesh_scalar)
print(mesh_scalar.array_names)

p = pv.Plotter()
p.add_mesh(mesh_scalar,scalars="u",cmap=colormap)
p.show_grid()
p.show()