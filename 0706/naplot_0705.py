import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 24})
data = np.loadtxt("na_0705.txt")

plt.clf()
plt.plot(data[:,1],data[:,2],lw=2,label=r"$|111\rangle\langle 111|$")

plt.plot(data[:,1],data[:,3],lw=2,label=r"$|r11\rangle\langle r11|$")
plt.plot(data[:,1],data[:,4],lw=2,label=r"$|1r1\rangle\langle 1r1|$")

plt.plot(data[:,1],data[:,5],lw=2,label=r"$|11r\rangle\langle 11r|$")
plt.plot(data[:,1],data[:,6],lw=2,label=r"$|rr1\rangle\langle rr1|$")

plt.plot(data[:,1],data[:,7],lw=2,label=r"$|1r1\rangle\langle 1r1|$")
plt.plot(data[:,1],data[:,8],lw=2,label=r"$|1rr\rangle\langle 1rr|$")
plt.plot(data[:,1],data[:,9],lw=2,label=r"$|rrr\rangle\langle rrr|$")


plt.xlabel("Time")
plt.xlim([0,6])
plt.yscale('log')
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig("0706_111.png",dpi=300)
plt.show()

