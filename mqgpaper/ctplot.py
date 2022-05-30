from tabulate import tabulate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 12})

dataarp = np.loadtxt("arp_3_ct.txt")
datasp = np.loadtxt("sp_3_F_ct.txt")

fig,ax=plt.subplots(1,1)
X=dataarp[:,3]
Y=dataarp[:,4]
Z=dataarp[:,2]
cp = ax.tricontourf(X, Y, Z)

fig.colorbar(cp)
ax.set_title(r"$C_kZ$ Gate Fidelity")
ax.set_xlabel('Control Qubit Phase Deviation')
ax.set_ylabel('1st Target Qubit Phase Deviation')
plt.savefig("arp_ct.pdf",dpi=300)
plt.show()
