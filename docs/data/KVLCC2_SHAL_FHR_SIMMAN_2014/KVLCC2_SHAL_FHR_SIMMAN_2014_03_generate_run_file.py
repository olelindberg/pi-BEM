import pandas as pd
numMeshes = 8
df_fs = pd.read_csv('pmm_app_shal_fhr_App01_pure_surge_full_scale.csv')

testNames = df_fs['testId'].values
numTests = len(testNames)
#mpiCommand = "/usr/bin/mpirun"
mpiCommand = "mpiexec"
file = open("runKVLCC2.sh", "w")
for testId in range(0, numTests):

    testName = testNames[testId]

    for meshId in range(0, numMeshes):

        meshName = "mesh0" + str(meshId)
        pathName = testName + "/" + meshName
        command = mpiCommand + " -np 4 ./bem_fma_3d 4 --input-path" + \
            " ../docs/data/KVLCC2_SHAL_FHR_SIMMAN_2014/" + pathName + "\n"
        print(command)
        file.write(command)
file.close()
