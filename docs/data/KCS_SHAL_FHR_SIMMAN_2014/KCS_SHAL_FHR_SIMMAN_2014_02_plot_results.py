import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

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

Tm_model = 0.2050
Tm_full = 10.8
df = pd.read_csv('pmm_app_shal_fhr_App01_pure_surge.csv')
df_fs = pd.read_csv('pmm_app_shal_fhr_App01_pure_surge_full_scale.csv')

U_fullscale = df_fs['surge velocity [m/s]'].values
h_fullscale = df_fs['depth [m]'].values
Frh_fullscale = U_fullscale/np.sqrt(g*h_fullscale)
FrL_fullscale = U_fullscale/np.sqrt(g*Lpp_full)
print(Frh_fullscale)
print(FrL_fullscale)


testIds = df['testId'].values
U = df['surge velocity [m/s]'].values
h = df['depth [m]'].values
sink_m = df['sink_m [mm]'].values
sink_a = df['sink_a [mm]'].values
sink_f = df['sink_f [mm]'].values
trim = df['trim [mm/m]'].values


Frh = U/np.sqrt(g*h)
FrL = U/np.sqrt(g*Lpp_model)
print(Frh)
print(FrL)

heaveforceAll = []
pitchmomentAll = []
for meshId in range(1, 9):

    heaveforce = []
    pitchmoment = []
    for testId in testIds[0:1]:

        filename = "mesh0" + str(meshId) + "/" + testId + "/output/force.csv"
        print(filename)
        data = np.genfromtxt(filename, delimiter=",")

        heaveforce.append(data[0][2])
        pitchmoment.append(data[0][4])

    heaveforceAll.append(heaveforce)
    pitchmomentAll.append(pitchmoment)

heaveforceAll = np.array(heaveforceAll)
pitchmomentAll = np.array(pitchmomentAll)

print(heaveforceAll)


if False:
    plt.figure()
    plt.plot(Frh[perm], 0.001*sink_m[perm]/Tm_model, 'bo-', label="EFD")
    plt.plot(Frh[perm], -np.array(sinkage[perm]) /
             Tm_full, 'ro-', label="BEM")
    plt.plot(Frh[perm], np.array(sinkage_barrass_restricted[perm]) /
             Tm_full, 'go-', label="Barrass, restricted, S=0.25")
    plt.plot(Frh[perm], np.array(sinkage_barrass_unrestricted[perm]) /
             Tm_full, 'co-', label="Barrass, unrestricted, S=0.1")
    plt.xlabel(r"$Fr_h []$")
    plt.ylabel(r"$sinkage/T_m [\%]$")
    plt.grid(True)
    plt.legend()

plt.figure(1)
plt.plot(range(0, len(heaveforceAll[:, 0])),
         heaveforceAll[:, 0], 'ro-', label="BEM")
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.grid(True)
plt.legend()


plt.figure(2)
plt.plot(range(0, len(pitchmomentAll[:, 0])),
         pitchmomentAll[:, 0], 'ro-', label="BEM")
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.grid(True)
plt.legend()


plt.show()
