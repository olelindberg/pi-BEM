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

quantile_fraction = 0.95

for i in range(34):

    # Read the source file.
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName("/home/ole/dev/projects/pi-BEM/docs/data/KCS_Limfjord/output/error_estimators" + str(i) + ".vtu")
    reader.Update()  # Needed because of GetScalarRange
    output = reader.GetOutput()

    #error_estimate = VN.vtk_to_numpy(output.GetPointData().GetArray('error_estimator_potential'))
    error_estimate = VN.vtk_to_numpy(output.GetPointData().GetArray('error_estimator_velocity'))

    x   = np.sort(error_estimate)
    mu  = np.sum(np.log(x))/len(x)
    s2  = np.sum((np.log(x)-mu)**2)/len(x)
    s   = np.sqrt(s2)   
    pdf = 1/(x*s*np.sqrt(2*np.pi))*np.exp(-((np.log(x)-mu)**2)/(2*s2))

    plt.figure()
    counts, bins = np.histogram(error_estimate,bins=2**6)
    plt.hist((bins[:-1]), (bins), weights=counts/np.sum(np.diff(bins)*counts))
    plt.axvline(x = np.exp(mu+s2/2), color = 'r', label = 'mean')
    plt.axvline(x = np.exp(mu+s2/2-3*s2), color = 'g', label = 'mean - x*var')
    plt.axvline(x = np.exp(mu+s2/2+3*s2), color = 'b', label = 'mean + x*var')
    plt.plot((x),pdf)
    plt.legend()
    plt.grid(True)    


plt.show() 


