import numpy as np 
import matplotlib.pyplot as plt

savefigures = False

g               = 9.80665
density         = 1000
Lpp_full        = 355
Lpp_model       = 3.984
Tm_model        = 0.163
Tm_full         = 14.5
waterplanearea  = 15315 # 22032
depth_model     = np.array([2.0*Tm_model,2.0*Tm_model,1.2*Tm_model])
velocity_model  = [0.327,0.872,0.327]
FnL              = velocity_model/np.sqrt(g*Lpp_model)
Fnh              = velocity_model/np.sqrt(g*depth_model)
velocity_full   = FnL*np.sqrt(g*Lpp_full)

F1 = np.genfromtxt("C1/output/force.csv",delimiter=",")
F2 = np.genfromtxt("C2/output/force.csv",delimiter=",")
F3 = np.genfromtxt("C3/output/force.csv",delimiter=",")
print(F1)
print(F2)
print(F3)
heaveforce = np.array([F1[-1],F2[-1],F3[0][2]])
print(heaveforce)
sinkage         = heaveforce/(g*density*waterplanearea)
print(sinkage)

C1 = np.genfromtxt("C1.dat",delimiter="")
C2 = np.genfromtxt("C2.dat",delimiter="")
C3 = np.genfromtxt("C3.dat",delimiter="")

t = np.array([0,230])

mm_to_m = 0.001
plt.figure(1)
plt.plot(C1[:,0],-C1[:,5]*mm_to_m/Tm_model*100,label="EFD",color='k')
plt.plot(t,0*t+sinkage[0]/Tm_full*100,label='BEM',color='r')
plt.title(r'C1: $Fn_L$=' + "{:1.2f}".format(FnL[0]) + r', $Fn_h$=' + "{:1.2f}".format(Fnh[0]))
plt.xlabel(r'$t~[s]$')
plt.ylabel(r'$heave/T_m ~ [\%]$')
plt.grid(True)
plt.legend()
plt.xlim([0,230])
if savefigures:
    filename = 'DTC_MASHCON_2019_C1.pdf'
    plt.savefig(filename,bbox_inches='tight')
    filename = 'DTC_MASHCON_2019_C1.png'
    plt.savefig(filename,bbox_inches='tight')

plt.figure(2)
plt.plot(C2[:,0],-C2[:,5]*mm_to_m/Tm_model*100,label="EFD",color='k')
plt.plot(t,0*t+sinkage[1]/Tm_full*100,label='BEM',color='r')
plt.title(r'C2: $Fn_L$=' + "{:1.2f}".format(FnL[1]) + r', $Fn_h$=' + "{:1.2f}".format(Fnh[1]))
plt.xlabel(r'$t~[s]$')
plt.ylabel(r'$heave/T_m ~ [\%]$')
plt.grid(True)
plt.legend()
plt.xlim([0,110])
if savefigures:
    filename = 'DTC_MASHCON_2019_C2.pdf'
    plt.savefig(filename,bbox_inches='tight')
    filename = 'DTC_MASHCON_2019_C2.png'
    plt.savefig(filename,bbox_inches='tight')

plt.figure(3)
plt.plot(C3[:,0],-C3[:,5]*mm_to_m/Tm_model*100,label="EFD",color='k')
plt.plot(t,0*t+sinkage[2]/Tm_full*100,label='BEM',color='r')
plt.title(r'C3: $Fn_L$=' + "{:1.2f}".format(FnL[2]) + r', $Fn_h$=' + "{:1.2f}".format(Fnh[2]))
plt.xlabel(r'$t~[s]$')
plt.ylabel(r'$heave/T_m ~ [\%]$')
plt.grid(True)
plt.legend()
plt.xlim([0,230])
if savefigures:
    filename = 'DTC_MASHCON_2019_C3.pdf'
    plt.savefig(filename,bbox_inches='tight')
    filename = 'DTC_MASHCON_2019_C3.png'
    plt.savefig(filename,bbox_inches='tight')

plt.show()

