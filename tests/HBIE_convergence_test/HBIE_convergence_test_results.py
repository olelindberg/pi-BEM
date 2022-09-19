


from cProfile import label
import matplotlib.pyplot as plt
import numpy as np

show_ndof_vs_err = True
show_ndof_vs_cputime = True
show_cputime_vs_error = True

#data = np.genfromtxt("tests/HBIE_convergence_test/HBIE_convergence_test1.csv",delimiter=",")
data = np.genfromtxt("/home/ole/dev/projects/pi-BEM/build_ole_hsing_coeffs/tests/HBIE_convergence_test/HBIE_convergence_test.csv",delimiter=",")

print(data)
ndof = data[:,1]
order = data[:,2]
err_l2 = data[:,-3]
err_inf = data[:,-2]
time = data[:,-1]
order_max = 4


R    = 1
area = 4*np.pi*R**2
element_area = area/ndof # actually not the element error, more like the dof area
dx = np.sqrt(element_area)

if (show_ndof_vs_err):
    dxx = np.sort(np.unique(dx))
    plt.plot(np.log10(dxx),1*np.log10(dxx),'k--')
    plt.plot(np.log10(dxx),2*np.log10(dxx),'k--')
    plt.plot(np.log10(dxx),3*np.log10(dxx),'k--')
    plt.plot(np.log10(dxx),4*np.log10(dxx),'k--')
    plt.plot(np.log10(dxx),5*np.log10(dxx),'k--')

    for p in range(1,order_max+1):
        id = np.where(order==p)
        lbl = "degree = " + str(p)
        plt.plot(np.log10(dx[id]),np.log10(err_inf[id]),'-o',label=lbl + ", Linf")
        plt.plot(np.log10(dx[id]),np.log10(err_l2[id]),'-o',label=lbl + ", L2")


    plt.grid(True)
    plt.xlabel(r"$\log_{10}(dx)$")
    plt.ylabel(r"$\log_{10}(||err||)$")
    plt.legend()

if (show_ndof_vs_cputime):
    plt.figure()
    for p in range(1,order_max+1):
        id = np.where(order==p)
        lbl = "degree = " + str(p)
        plt.plot(np.log10(ndof[id]),np.log10(time[id]),'-o',label=lbl)

    plt.grid(True)
    plt.xlabel(r"$\log_{10}(ndof)$")
    plt.ylabel(r"$\log_{10}(cputime[s])$")
    plt.legend()


if (show_cputime_vs_error):
    plt.figure()
    for p in range(1,order_max+1):
        id = np.where(order==p)
        lbl = "degree = " + str(p)
        plt.plot(np.log10(time[id]),np.log10(err_inf[id]),'-o',label=lbl + ", Linf")
        plt.plot(np.log10(time[id]),np.log10(err_l2[id]),'-o',label=lbl + ", L2")

    plt.grid(True)
    plt.xlabel(r"$\log_{10}(cputime [s])$")
    plt.ylabel(r"$\log_{10}(||err||)$")
    plt.legend()


plt.show()