import matplotlib.pyplot as plt
import numpy as np
import math

savefigures         = True

colors = 'rgbcmyk'
grav   = 9.80665
LPP    = 2.5
Tm     = 0.156
depth  = 0.24
Fn_L   = np.array([0.05, 0.08, 0.11, 0.15, 0.18, 0.21, 0.24])
vel    = Fn_L*math.sqrt(grav*LPP)
Fn_h   = vel/math.sqrt(grav*depth)

heaveAll  = np.zeros((7,1))
velxAll  = np.zeros((7,1))

cfddata    = np.genfromtxt('WigleyHull_CFDResults.csv',delimiter=',')
cfdFr      = cfddata[:,1]
cfdsinkage = cfddata[:,4]

efddata    = np.genfromtxt('WigleyHullSinkageFig10.csv',delimiter=',')
efdFr      = efddata[:,1]
efdsinkage = efddata[:,2]


bswdata    = np.genfromtxt('WigleyHullSinkageBoussinesq.csv',delimiter=',')
bswFr      = bswdata[:,1]
bswsinkage = bswdata[:,2]

bemdata    = np.genfromtxt('WigleyHullSinkage_piBEM.csv',delimiter=',')
bemFr      = bemdata[:,1]
bemsinkage = bemdata[:,2]

plt.figure(5)
plt.plot(Fn_h, efdsinkage/Tm*100.0, 'k-o',label='EFD')
plt.plot(Fn_h,  bswsinkage/Tm*100.0, 'r-^',label='BSW')
plt.plot(Fn_h,  bemsinkage/Tm*100.0, 'g-^',label='BEM')
plt.plot(Fn_h,-cfdsinkage/Tm*100.0, 'b-s',label='CFD')
plt.grid(True)
plt.legend()
plt.xlabel(r'$Fn_h$' + ' ' + r'$[]$')
plt.ylabel(r'$sinkage/T_m$' + ' ' + r'$[\%]$')
plt.xlim(0.1,0.8)
plt.ylim(0,7.0)

if savefigures:
    filename = 'Beguin_MASHCON2013_WigleyHull_Sinkage.pdf'
    plt.savefig(filename,bbox_inches='tight')
    filename = 'Beguin_MASHCON2013_WigleyHull_Sinkage.png'
    plt.savefig(filename,bbox_inches='tight')

plt.show()

