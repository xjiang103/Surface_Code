import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 24})
data = np.loadtxt("na_1108_1.txt")

plt.clf()
plt.plot(data[:,1],data[:,2],lw=2,label=r"$|1111\rangle\langle 1111|$")

plt.plot(data[:,1],data[:,3],lw=2,label=r"$|r111\rangle\langle r111|$")
plt.plot(data[:,1],data[:,4],lw=2,label=r"$|1r11\rangle\langle 1r11|$")

plt.plot(data[:,1],data[:,5],lw=2,label=r"$|rr11\rangle\langle rr11|$")
plt.plot(data[:,1],data[:,6],lw=2,label=r"$|1rr1\rangle\langle 1rr1|$")

plt.plot(data[:,1],data[:,7],lw=2,label=r"$|rrr1\rangle\langle rrr1|$")
plt.plot(data[:,1],data[:,8],lw=2,label=r"$|1rrr\rangle\langle 1rrr|$")

#plt.plot(data[:,1],data[:,8]+3*data[:,7]+3*data[:,6]+3*data[:,5]+3*data[:,4]+data[:,3],lw=2,label="sum")

plt.xlabel("Time")
plt.xlim([0,0.55])
#plt.yscale('log')
plt.legend(loc='upper right')
plt.tight_layout()
plt.show()
plt.savefig("pops.png",dpi=300)
