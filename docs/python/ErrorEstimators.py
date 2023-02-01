from vtk.util import numpy_support as VN
import vtkmodules.all as vtk
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from fitter import Fitter, get_common_distributions, get_distributions


colors = vtk.vtkNamedColors()

for i in range(38):

    # Read the source file.
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName("/home/ole/dev/projects/pi-BEM/docs/data/KCS_Limfjord/output/error_estimators" + str(i) + ".vtu")
    reader.Update()  # Needed because of GetScalarRange
    output = reader.GetOutput()

    u = VN.vtk_to_numpy(output.GetPointData().GetArray('error_estimator_potential'))
    b = VN.vtk_to_numpy(output.GetPointData().GetArray('error_estimator_velocity'))

    f = Fitter(u,distributions= get_distributions())
    f.fit()
    print("potential")
    print(f.get_best(method = 'sumsquare_error'))
    plt.figure()
    f.summary()

    f = Fitter(b,distributions= get_distributions())
    f.fit()
    print("velocity")
    print(f.get_best(method = 'sumsquare_error'))
    plt.figure()
    f.summary()

plt.show() 


