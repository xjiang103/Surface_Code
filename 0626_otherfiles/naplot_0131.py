import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 24})
data = np.loadtxt("na_0131.txt")

plt.clf()
plt.plot(data[:,1],data[:,2],lw=2,label=r"$|1111\rangle\langle 1111|$")
##
plt.plot(data[:,1],data[:,3],lw=2,label=r"$|1110\rangle\langle 1110|$")
plt.plot(data[:,1],data[:,4],lw=2,label=r"$|1100\rangle\langle 1100|$")
##
##plt.plot(data[:,1],3*data[:,5],lw=2,label=r"$3*|rr11\rangle\langle rr11|$")
##plt.plot(data[:,1],3*data[:,6],lw=2,label=r"$3*|1rr1\rangle\langle 1rr1|$")
##
##plt.plot(data[:,1],3*data[:,7],lw=2,label=r"$3*|rrr1\rangle\langle rrr1|$")
##plt.plot(data[:,1],data[:,8],lw=2,label=r"$|1rrr\rangle\langle 1rrr|$")
##
##plt.plot(data[:,1],data[:,8]+3*data[:,7]+3*data[:,6]+3*data[:,5]+3*data[:,4]+data[:,3],lw=2,label="sum")

##plt.plot(data[:,1],data[:,9],lw=2,label="phase of state 1111")
####plt.plot(data[:,1],data[:,10],lw=2,label="phase of state 1110")
####plt.plot(data[:,1],data[:,11],lw=2,label="phase2")
####plt.plot(data[:,1],data[:,12],lw=2,label="phase3")
####plt.plot(data[:,1],data[:,13],lw=2,label="phase4")
##plt.axhline(y = 0, color ="green", linestyle ="--")
##plt.axhline(y = -0.2021, color ="red", linestyle =":",label="-0.2021")
##plt.axhline(y = 3.14, color ="green", linestyle ="--")
##plt.axhline(y = -3.14, color ="green", linestyle ="--")

plt.xlabel("Time")
plt.xlim([0,2])
#plt.yscale('log')
plt.legend(loc='upper right')
plt.tight_layout()
plt.show()
plt.savefig("pops.png",dpi=300)
