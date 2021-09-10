import numpy as np 
import matplotlib.pyplot as plt


g               = 9.80665
density         = 1000
Lpp_full        = 355
Lpp_model       = 3.984
scale           = Lpp_full/Lpp_model
Tm_model        = 0.163
Tm_full         = 14.5
waterplanearea  = 22032
depth_model     = np.array([2.0*Tm_model,2.0*Tm_model,1.2*Tm_model])
velocity_model  = [0.327,0.872,0.327]
FnL              = velocity_model/np.sqrt(g*Lpp_model)
Fnh              = velocity_model/np.sqrt(g*depth_model)
velocity_full   = FnL*np.sqrt(g*Lpp_full)

C1 = np.genfromtxt("C1.dat",delimiter="")
C2 = np.genfromtxt("C2.dat",delimiter="")
C3 = np.genfromtxt("C3.dat",delimiter="")

t = np.array([0,230])

mm_to_m = 0.001
sampling_freq = 0.05

i1 = int(125/sampling_freq)
i2 = int(175/sampling_freq)
zFPP_avg = np.average(C1[i1:i2,8])
zAPP_avg = np.average(C1[i1:i2,9])

print("C1 zAPP_avg full " + str(zAPP_avg))
print("C1 zFPP_avg full " + str(zFPP_avg))
dz       = zFPP_avg-zAPP_avg
trim     = dz/Lpp_model
print("C1 trim " + str(np.rad2deg(np.arcsin(trim*mm_to_m))) + " deg")

plt.figure(1)
plt.plot(C1[:,0],-C1[:,5]*mm_to_m,label="heave",color='r')
plt.plot(C1[:,0],-C1[:,8]*mm_to_m,label="z FPP",color='g')
plt.plot(C1[:,0],-C1[:,9]*mm_to_m,label="z APP",color='b')
plt.plot(C1[:,0],0*C1[:,0]-zFPP_avg*mm_to_m,label="z FPP avg",color='m')
plt.plot(C1[:,0],0*C1[:,0]-zAPP_avg*mm_to_m,label="z APP avg",color='c')
plt.title(r'C1: $Fn_L$=' + "{:1.2f}".format(FnL[0]) + r', $Fn_h$=' + "{:1.2f}".format(Fnh[0]))
plt.xlabel(r'$t~[s]$')
plt.ylabel(r'$heave/T_m ~ [\%]$')
plt.grid(True)
plt.legend()
plt.xlim([0,230])

plt.figure(2)
plt.plot(C1[:,0],C1[:,6],label="trim",color='r')
plt.plot(C1[:,0],0*C1[:,0]-trim,label="trim avg",color='m')
plt.title(r'C1: $Fn_L$=' + "{:1.2f}".format(FnL[0]) + r', $Fn_h$=' + "{:1.2f}".format(Fnh[0]))
plt.xlabel(r'$t~[s]$')
plt.ylabel(r'$trim ~ [mm/m]$')
plt.grid(True)
plt.legend()
plt.xlim([0,230])



i1 = int(40/sampling_freq)
i2 = int(70/sampling_freq)
zFPP_avg = np.average(C2[i1:i2,8])
zAPP_avg = np.average(C2[i1:i2,9])
print("C2 zAPP_avg full " + str(zAPP_avg))
print("C2 zFPP_avg full " + str(zFPP_avg))
dz       = zFPP_avg-zAPP_avg
trim     = dz/Lpp_model
print("C2 trim " + str(np.rad2deg(np.arcsin(trim*mm_to_m))) + " deg")


plt.figure(3)
plt.plot(C2[:,0],-C2[:,5]*mm_to_m,label="heave",color='r')
plt.plot(C2[:,0],-C2[:,8]*mm_to_m,label="z FPP",color='g')
plt.plot(C2[:,0],-C2[:,9]*mm_to_m,label="z APP",color='b')
plt.plot(C2[:,0],0*C2[:,0]-zFPP_avg*mm_to_m,label="z FPP avg",color='m')
plt.plot(C2[:,0],0*C2[:,0]-zAPP_avg*mm_to_m,label="z APP avg",color='c')
plt.title(r'C2: $Fn_L$=' + "{:1.2f}".format(FnL[1]) + r', $Fn_h$=' + "{:1.2f}".format(Fnh[1]))
plt.xlabel(r'$t~[s]$')
plt.ylabel(r'$heave/T_m ~ [\%]$')
plt.grid(True)
plt.legend()
plt.xlim([0,110])

plt.figure(4)
plt.plot(C2[:,0],C2[:,6],label="trim",color='r')
plt.plot(C2[:,0],0*C2[:,0]-trim,label="trim avg",color='m')
plt.title(r'C2: $Fn_L$=' + "{:1.2f}".format(FnL[1]) + r', $Fn_h$=' + "{:1.2f}".format(Fnh[1]))
plt.xlabel(r'$t~[s]$')
plt.ylabel(r'$trim ~ [mm/m]$')
plt.grid(True)
plt.legend()
plt.xlim([0,230])


i1 = int(125/sampling_freq)
i2 = int(175/sampling_freq)
zFPP_avg = np.average(C3[i1:i2,8])
zAPP_avg = np.average(C3[i1:i2,9])
print("C3 zAPP_avg full " + str(scale*zAPP_avg))
print("C3 zFPP_avg full " + str(scale*zFPP_avg))
dz       = zFPP_avg-zAPP_avg
trim     = dz/Lpp_model
print("C3 trim " + str(np.rad2deg(np.arcsin(trim*mm_to_m))) + " deg")

plt.figure(5)
plt.plot(C3[:,0],-C3[:,5]*mm_to_m,label="heave",color='r')
plt.plot(C3[:,0],-C3[:,8]*mm_to_m,label="z FPP",color='g')
plt.plot(C3[:,0],-C3[:,9]*mm_to_m,label="z APP",color='b')
plt.plot(C3[:,0],0*C3[:,0]-zFPP_avg*mm_to_m,label="z FPP avg",color='m')
plt.plot(C3[:,0],0*C3[:,0]-zAPP_avg*mm_to_m,label="z APP avg",color='c')
plt.title(r'C3: $Fn_L$=' + "{:1.2f}".format(FnL[2]) + r', $Fn_h$=' + "{:1.2f}".format(Fnh[2]))
plt.xlabel(r'$t~[s]$')
plt.ylabel(r'$heave/T_m ~ [\%]$')
plt.grid(True)
plt.legend()
plt.xlim([0,230])

plt.figure(6)
plt.plot(C3[:,0],C3[:,6],label="trim",color='r')
plt.plot(C3[:,0],0*C3[:,0]-trim,label="trim avg",color='m')
plt.title(r'C3: $Fn_L$=' + "{:1.2f}".format(FnL[2]) + r', $Fn_h$=' + "{:1.2f}".format(Fnh[2]))
plt.xlabel(r'$t~[s]$')
plt.ylabel(r'$trim ~ [mm/m]$')
plt.grid(True)
plt.legend()
plt.xlim([0,230])

plt.show()
