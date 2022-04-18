import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 24})
data = np.loadtxt("phase_0105.txt")

plt.clf()
plt.plot(data[:,0],data[:,1],lw=2,label=r"$|1111\rangle\langle 1111|$")

#plt.plot(data[:,0],data[:,2],lw=2,label=r"$|1110\rangle\langle 1110|$")

#plt.plot(data[:,0],data[:,1]-data[:,2],lw=2,label=r"$|1111\rangle\langle 1111|$-$|1110\rangle\langle 1110|$")

plt.axhline(y = 0, color ="green", linestyle ="--")
plt.axhline(y = -0.2021, color ="red", linestyle =":",label="-0.2021")
plt.axhline(y = 3.14, color ="green", linestyle ="--")
plt.axhline(y = -3.14, color ="green", linestyle ="--")

plt.xlabel("Time")
plt.xlim([0,2])
#plt.yscale('log')
plt.legend(loc='upper right')
plt.tight_layout()
plt.show()
plt.savefig("pops.png",dpi=300)
