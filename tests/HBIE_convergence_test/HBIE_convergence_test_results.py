from cProfile import label
import matplotlib.pyplot as plt
import numpy as np

data = np.genfromtxt("/home/ole/Projects/pi-BEM/tests/HBIE_convergence_test/HBIE_convergence_test_results.dat",delimiter=",")
print(data)

R    = 1
area = 4*np.pi*R**2
element_area = area/data[:,4] # actually not the element error, more like the dof area
element_length = np.sqrt(element_area)

plt.plot(np.log10(element_length[:4]),0.5+1*np.log10(element_length[:4]),'k--',label="slope 1")
plt.plot(np.log10(element_length),np.log10(data[:,5]),'ro',label='$L_{\infty} error$')
plt.plot(np.log10(element_length),np.log10(data[:,6]),'go',label='$L_2$ error')
plt.grid(True)
plt.xlabel(r"$\log_{10}(\Delta x)$")
plt.ylabel(r"$\log_{10}(||err||)$")
plt.legend()
plt.show()