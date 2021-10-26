import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

showsinkage = True

g = 9.80665
density = 1000
scale = 75
Lpp_full = 320
Tm_full = 20.8
waterplanearea = 16963.3
V0 = 312600
zG = 20.8

Lpp_model = 1/scale*Lpp_full
Tm_model = 1/scale*Tm_full


df = pd.read_csv('pmm_app_shal_fhr_App01_pure_surge.csv')
df_fs = pd.read_csv('pmm_app_shal_fhr_App01_pure_surge_full_scale.csv')

U_fullscale = df_fs['surge velocity [m/s]'].values
h = df_fs['depth [m]'].values
Frh = U_fullscale/np.sqrt(g*h)
FrL = U_fullscale/np.sqrt(g*Lpp_full)
print(Frh)
print(FrL)


testIds = df['testId'].values
U = df['surge velocity [m/s]'].values
h = df['depth [m]'].values
sink_m = df['sink_m [mm]'].values

Frh = U/np.sqrt(g*h)
FrL = U/np.sqrt(g*Lpp_model)
print(Frh)
print(FrL)

perm = np.argsort(Frh)
if (showsinkage):
    plt.figure(1)
    plt.plot(Frh[perm], 100*0.001*sink_m[perm] /
             Tm_model, 'ko-', label='EFD')


for gridId in range(1, 4):

    sinkage = []
    pitchmoment = []

    for testId in testIds:
        filename = "pibem/grid0" + \
            str(gridId) + "/" + testId + "/output/force.csv"
        print(filename)
        data = np.genfromtxt(filename, delimiter=",")
        heaveforce = data[0][2]
        sinkage.append(heaveforce/(g*density*waterplanearea))

        pitchmoment.append(data[0][4])

    sinkage = np.array(sinkage)
    pitchmoment = np.array(pitchmoment)

    #pitch = pitchmoment/((g*density*V0*(zB-zG+Sxx/V0)))

    if showsinkage:
        plt.figure(1)
        plt.plot(Frh[perm], -100*np.array(sinkage[perm]) /
                 Tm_full, 'o-', label="grid " + str(gridId))
        plt.xlabel(r"$Fr_h []$")
        plt.ylabel(r"$sinkage/T_m [\%]$")
        plt.grid(True)
        plt.legend()

    plt.figure(2)
  #  plt.plot(Frh[perm], 0.001*trim[perm], 'bo-', label="EFD")
    plt.plot(Frh[perm], np.array(pitchmoment[perm]),
             'o-', label="grid " + str(gridId))
    plt.xlabel(r"$Fr_h~ []$")
    plt.ylabel(r"$pitch~moment~[Nm]$")
    plt.grid(True)
    plt.legend()


plt.show()
