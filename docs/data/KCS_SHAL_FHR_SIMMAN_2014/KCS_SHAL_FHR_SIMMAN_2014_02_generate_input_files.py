import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import shutil
from string import Template

g = 9.80665
density = 1000
scale = 52.667
Lpp_full = 230
Lpp_model = 1/scale*Lpp_full
CB = 0.651
zB = 5.91163
zG = 11.4
Sxx = 467132
Syy = 1.96861e+07
A = 6227.87
V0 = 52030

Tm = 10.8
df_fs = pd.read_csv('pmm_app_shal_fhr_App01_pure_surge_full_scale.csv')

velocities = df_fs['surge velocity [m/s]'].values
testNames = df_fs['testId'].values
UKCs = df_fs['depth under keel [m]'].values
numTests = len(UKCs)
print(testNames)
print(UKCs)
print(velocities)

for meshId in range(1, 9):
    meshName = "mesh0" + str(meshId)
    for testId in range(0, numTests):

        testName = testNames[testId]
        vel = velocities[testId]
        ukc = UKCs[testId]
        pathName = meshName + "/" + testName
        print(pathName)
        print(vel)
        print(ukc)
        shutil.copytree("pibem_setup_template", pathName)

        # Parameter file:
        filename = pathName + "/parameters_bem_3.prm"
        file = open(filename, "r")
        text = file.read()
        file.close()

        temp_obj = Template(text)
        dct = {'velocity': str(vel)}
        temp_obj = temp_obj.substitute(dct)
        text = str(temp_obj)

        file = open(filename, "w")
        file.write(text)
        file.close()

        # Mesh file:
        filename = pathName + "/mesh.inp"
        file = open(filename, "r")
        text = file.read()
        file.close()

        temp_obj = Template(text)
        dct = {'UKC': str(ukc)}
        temp_obj = temp_obj.substitute(dct)
        text = str(temp_obj)

        file = open(filename, "w")
        file.write(text)
        file.close()


file = open("runKCS.sh", "w")
for meshId in range(1, 9):
    meshName = "mesh0" + str(meshId)
    for testId in range(0, numTests):

        testName = testNames[testId]
        pathName = meshName + "/" + testName
        command = "mpiexec -np 4 ./bem_fma_3d 4 --input-path" + \
            " ../docs/data/KCS_SHAL_FHR_SIMMAN_2014/" + pathName + "\n"
        file.write(command)
file.close()
