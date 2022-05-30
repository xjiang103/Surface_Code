from tabulate import tabulate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 12})

dataarp = np.loadtxt("arp_3_ct.txt")
datasp = np.loadtxt("sp_3_F_ct.txt")

fig,ax=plt.subplots(1,1)
X=datasp[:,3]
Y=datasp[:,4]
Z=datasp[:,2]
cp = ax.tricontourf(X, Y, Z)

fig.colorbar(cp)
ax.set_title(r"$CZ_k$ Gate Fidelity")
ax.set_xlabel('Control Qubit Phase Deviation')
ax.set_ylabel('1st Target Qubit Phase Deviation')
plt.savefig("sp_ct.pdf",dpi=300)
plt.show()
