import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import shutil
from string import Template


def substitudeTemplate(filename, dct):

    file = open(filename, "r")
    text = file.read()
    file.close()

    temp_obj = Template(text)
    temp_obj = temp_obj.substitute(dct)
    text = str(temp_obj)

    file = open(filename, "w")
    file.write(text)
    file.close()
# ------------------------------------------------------------------------------


substitude = True
numMeshes = 8
df_fs = pd.read_csv('pmm_app_shal_fhr_App01_pure_surge_full_scale.csv')

testNames = df_fs['testId'].values
velocities = df_fs['surge velocity [m/s]'].values
UKCs = df_fs['depth under keel [m]'].values

numTests = len(testNames)
print(testNames)
print(UKCs)
print(velocities)
for testId in range(0, numTests):

    testName = testNames[testId]

    for meshId in range(0, numMeshes):

        meshName = "mesh0" + str(meshId)

        pathName = testName + "/" + meshName
        vel = velocities[testId]
        ukc = UKCs[testId]
        cellSizeMin = 1/1.51**(meshId)

        print("pathname      : " + pathName)
        print("velocity      : " + str(vel))
        print("UKC           : " + str(ukc))
        print("cell size min : " + str(cellSizeMin))

        if substitude:

            shutil.copytree("template", pathName)

            # Parameter file:
            filename = pathName + "/parameters_bem_3.prm"
            dct = {"velocity": str(vel)}
            substitudeTemplate(filename, dct)

            # Mesh file:
            filename = pathName + "/mesh.inp"
            dct = {'UKC': str(ukc)}
            substitudeTemplate(filename, dct)

            # pibem.json file:
            filename = pathName + "/pibem.json"
            dct = {'cellSizeMin': str(cellSizeMin)}
            substitudeTemplate(filename, dct)
