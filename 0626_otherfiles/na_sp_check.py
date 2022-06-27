from qutip import *
from qutip.qip.operations import cnot,snot

import numpy as np

#order of enumeration: 0 d 1 r
zero=basis(4,0)
one =basis(4,2)


cz_arp = np.diag([1,-1,-1,-1])
psi0 = tensor(basis(2,1),basis(2,1))
psi0 = tensor(snot(),snot())*psi0
psi0 = cz_arp*psi0

psi10 = ket2dm(Qobj(psi0))
psi10=Qobj(psi10,dims=[[2,2],[2,2]])
print(psi10)

psi0 = tensor(qeye(2),snot())*psi0
psi0 = ket2dm(Qobj(psi0))
psi0=Qobj(psi0,dims=[[2,2],[2,2]])

print(psi0)

mat_1=np.array([[0.25,0.1779-1j*0.1756,0.1779-1j*0.1756,-0.0036+1j*0.25],
                 [0.1779+1j*0.1756,0.25,0.25,-0.1782+1j*0.1754],
                 [0.1779+1j*0.1756,0.25,0.25,-0.1782+1j*0.1754],
                 [-0.0036-1j*0.25,-0.1782-1j*0.1754,-0.1782-1j*0.1754,0.25]])
d_r=Qobj(mat_1,dims=[[2,2],[2,2]])

print(d_r)

print(fidelity(psi10,d_r))
