import numpy as np 
import matplotlib.pyplot as plt

g               = 9.80665
density         = 1000

Lpp_full        = 355
Lpp_model       = 3.984
Tm_model        = 0.163
Tm_full         = 14.5
waterplanearea  = 22032

velocity_model  = [0.327,0.872,0.327]
Fn              = velocity_model/np.sqrt(g*Lpp_model)
velocity_full   = Fn*np.sqrt(g*Lpp_full)
heaveforce      = np.array([-3.39038e+06,-2.41094e+07])
sinkage         = heaveforce/(g*density*waterplanearea)
print(sinkage)
print(Fn)
print(velocity_full)

C1 = np.genfromtxt("C1.dat",delimiter="")
C2 = np.genfromtxt("C2.dat",delimiter="")
C3 = np.genfromtxt("C3.dat",delimiter="")

t = np.array([0,230])

mm_to_m = 0.001
plt.plot(C1[:,0],-C1[:,5]*mm_to_m/Tm_model)
plt.plot(C2[:,0],-C2[:,5]*mm_to_m/Tm_model)
plt.plot(C3[:,0],-C3[:,5]*mm_to_m/Tm_model)
plt.plot(t,0*t+sinkage[0]/Tm_full)
plt.plot(t,0*t+sinkage[1]/Tm_full)

plt.grid(True)
plt.xlim([0,230])
plt.show()

