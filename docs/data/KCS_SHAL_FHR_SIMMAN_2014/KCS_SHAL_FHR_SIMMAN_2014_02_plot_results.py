import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd

g = 9.80665
density = 1000
scale = 52.667
Lpp_full  = 230
Lpp_model = 1/scale*Lpp_full


Tm_model = 0.2050
Tm_full = 10.8
df = pd.read_csv('pmm_app_shal_fhr_App01_pure_surge.csv')
df_fs = pd.read_csv('pmm_app_shal_fhr_App01_pure_surge_full_scale.csv')

U       = df_fs['surge velocity [m/s]'].values
h       = df_fs['depth [m]'].values
Frh = U/np.sqrt(g*h)
FrL = U/np.sqrt(g*Lpp_full)
print(Frh)
print(FrL)


testIds = df['testId'].values
U       = df['surge velocity [m/s]'].values
h       = df['depth [m]'].values
sink_m  = df['sink_m [mm]'].values

Frh = U/np.sqrt(g*h)
FrL = U/np.sqrt(g*Lpp_model)
print(Frh)
print(FrL)

A = 6227.87
sinkage = []
for testId in testIds:
    filename = testId + "/output/force.csv"
    print(filename)
    data        = np.genfromtxt(filename,delimiter=",")
    heaveforce  = data[0][2]
    sinkage.append(heaveforce/(g*density*A))

sinkage = np.array(sinkage)

perm = np.argsort(Frh)
plt.figure()
plt.plot(Frh[perm],0.001*sink_m[perm]/Tm_model,'bo-')
plt.plot(Frh[perm],-np.array(sinkage[perm])/Tm_full,'ro-')
plt.xlabel(r"$Fr_h []$")
plt.ylabel(r"$sinkage/T_m [\%]$")
plt.grid(True)
plt.show()