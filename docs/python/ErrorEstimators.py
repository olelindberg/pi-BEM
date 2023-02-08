from distutils.log import error
from vtk.util import numpy_support as VN
import vtkmodules.all as vtk
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from fitter import Fitter, get_common_distributions, get_distributions
import math 

colors = vtk.vtkNamedColors()

for i in range(1):

    # Read the source file.
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName("/home/ole/dev/projects/pi-BEM/docs/data/KCS_Limfjord/output/errors" + str(i) + ".vtu")
    reader.Update()  # Needed because of GetScalarRange
    output = reader.GetOutput()

    error_estimate_pot = VN.vtk_to_numpy(output.GetPointData().GetArray('error_estimator_pot'))

    plt.figure()
    plt.plot(error_estimate_pot)
    plt.grid(True)    

    plt.figure()
    counts, bins = np.histogram(error_estimate_pot,bins=2**6)
    plt.hist((bins[:-1]), (bins), weights=counts/np.sum(np.diff(bins)*counts))
    plt.legend()
    plt.grid(True)    

    error_estimate_vel = VN.vtk_to_numpy(output.GetPointData().GetArray('error_estimator_vel'))

    plt.figure()
    plt.plot(error_estimate_vel)
    plt.grid(True)    

    plt.figure()
    counts, bins = np.histogram(error_estimate_vel,bins=2**6)
    plt.hist((bins[:-1]), (bins), weights=counts/np.sum(np.diff(bins)*counts))
    plt.legend()
    plt.grid(True)   

plt.show() 


