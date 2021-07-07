from qutip import *
from qutip.qip.operations import cnot,snot
from qutip.qip.circuit import QubitCircuit
from qutip.qip.gates import gate_sequence_product
import numpy as np

def user_gate1():
    diag=np.array([1,-1,-1,-1])
    mat=np.diag(diag)
    return Qobj(mat,dims=[[2,2],[2,2]])
def user_gate2():
    diag=np.array([-1,-1,-1,1])
    mat=np.diag(diag)
    return Qobj(mat,dims=[[2,2],[2,2]])                
cz_arp=user_gate2()

zero=basis(2,0)
one =basis(2,1)

data=rand_ket(4)
data=Qobj(data,dims=[[2,2],[1,1]])
#data=tensor(one,one)
data2=data

q=QubitCircuit(2,reverse_states=False)
q.add_gate("CNOT",controls=[0],targets=[1])
qmat=gate_sequence_product(q.propagators())

data=qmat*data

data2=tensor(qeye(2),snot())*data2
data2=cz_arp*data2
data2=tensor(qeye(2),snot())*data2

print(data)
print(data2)

print(fidelity(data,data2))
