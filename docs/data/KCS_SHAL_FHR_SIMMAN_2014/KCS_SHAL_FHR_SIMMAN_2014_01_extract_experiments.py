import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

g = 9.80665

sheetname = 'Environment'

df = pd.read_excel(
    'pmm_app_shal_fhr/test results/App01_KCS_PMM_results_May2014.xlsx', sheet_name=sheetname)
scale = df['Environment'][8]
LPP_full = df['Environment'][10]
LPP_model = df['Environment'][10]/scale
depth = df['hmodel (m)'][0:3].values
print("LPP_model    = " + str(LPP_model))
print("depth  = " + str(depth))
sheetnames = ['C0101A01', 'C0101A02', 'C0101A03']
cnt = 0

tests = []
for sheetname in sheetnames:

    df = pd.read_excel(
        'pmm_app_shal_fhr/test results/App01_KCS_PMM_results_May2014.xlsx', sheet_name=sheetname)
    df = df.fillna('')
    testname = sheetname + "_C"
    df = df.loc[df[df.columns[0]].str.contains(testname)]
    df = df.loc[df[df.columns[1]].str.contains("STATX0")]
    df = df[df[df.columns[4]] == 0]  # Zero drift angle
    df = df[df[df.columns[22]] == 0]  # Zero rudder angle

    testId = df[df.columns[0]].values
    Uc = df['VELOCITY'].values  # [mm]
    sink_m = df['DOF 3'].values  # [mm]
    trim = df['DOF 5'].values  # [mm/m]
    sink_a = sink_m-trim*LPP_model/2
    sink_f = sink_m+trim*LPP_model/2

    print('testId     = ' + str(testId))
    print('Uc     = ' + str(Uc))
    print('sink_m = ' + str(sink_m))

    for i in range(0, len(Uc)):
        tests.append([testId[i][0:11], depth[cnt], Uc[i],
                     sink_m[i], trim[i], sink_a[i], sink_f[i]])
    cnt = cnt+1


df = pd.DataFrame(tests)
df.columns = ['testId', 'depth [m]', 'surge velocity [m/s]',
              'sink_m [mm]', 'trim [mm/m]', 'sink_a [mm]', 'sink_f [mm]']
df.index = df['testId']
df = df.drop('testId', axis=1)
df.to_csv('pmm_app_shal_fhr_App01_pure_surge.csv')

Uc = df['surge velocity [m/s]'].values
FnL = Uc/np.sqrt(g*LPP_model)
print(FnL)

df['depth [m]'] = df['depth [m]']*scale
df['sink_m [mm]'] = df['sink_m [mm]']*scale
df['sink_a [mm]'] = df['sink_a [mm]']*scale
df['sink_f [mm]'] = df['sink_f [mm]']*scale
df['surge velocity [m/s]'] = FnL*np.sqrt(g*LPP_full)
df['depth under keel [m]'] = df['depth [m]']-10.8
df.to_csv('pmm_app_shal_fhr_App01_pure_surge_full_scale.csv')

Uc = df['surge velocity [m/s]'].values
FnL = Uc/np.sqrt(g*LPP_full)
print(FnL)
