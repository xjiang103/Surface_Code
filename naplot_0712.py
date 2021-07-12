import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
data = np.loadtxt("na_4_par_3lvl2_0712.txt")

plt.clf()
plt.plot(data[:,1],data[:,2],lw=2,label=r"$|1111\rangle\langle 1111|$")
plt.plot(data[:,1],data[:,3],lw=2,label=r"$|111r\rangle\langle 111r|$")
plt.plot(data[:,1],data[:,4],lw=2,label=r"$|11r1\rangle\langle 11r1|$")
plt.plot(data[:,1],data[:,5],lw=2,label=r"$|1r11\rangle\langle 1r11|$")
plt.plot(data[:,1],data[:,6],lw=2,label=r"$|r111\rangle\langle r111|$")

plt.xlabel("Time")
plt.xlim([0,0.54])
plt.yscale('log')
plt.legend(loc='upper right')
plt.tight_layout()
plt.show()
plt.savefig("pops.png",dpi=300)
