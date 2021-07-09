import pyvista as pv
import numpy as np
import math

#outputdir  ="example04a_hemisphere/output/"
outputdir  ="example04d_hemisphere_in_box_hT_1.2/output/"
#outputdir = "example04b_hemisphere_hT_2.0/output/"

filename    = outputdir + "result_scalar_results.vtu"
mesh_scalar = pv.read(filename)


filename    = "example04a_hemisphere/output/result_vector_results.vtu"
#filename    = "example04b_hemisphere_hT_2.0/output/result_vector_results.vtu"
filename    = outputdir + "result_vector_results.vtu"
mesh_vector = pv.read(filename)

pot       = mesh_scalar.get_array("phi")
vel       = mesh_vector.get_array(mesh_vector.array_names[0])
pnt       = mesh_vector.points

checkpoints = np.array([[0,0,-1],[0,-1,0],[0,1,0]])
checkpoints = np.array([[0,0,-1]])
checkpoints = np.array([[-1,0,0],[0,0,-1],[1,0,0]])
tol = 1e-3

print(checkpoints)
for j in range(0,checkpoints.shape[0]):
    x0 = checkpoints[j,0] 
    y0 = checkpoints[j,1]
    z0 = checkpoints[j,2]
    for i in range(0,vel.shape[0]):
        p = pnt[i]
        v = vel[i]
        pott = pot[i]

        if (abs(p[0]-x0)<1e-1 and abs(p[1]-y0)<tol and abs(p[2]-z0)<tol):
            print("\n")
            print(p)
            print(v)
 #           print(pott)            