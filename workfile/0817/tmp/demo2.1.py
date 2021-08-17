from IPython.display import Image
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
from qutip.qip.operations import cnot
from qutip.qip.circuit import QubitCircuit

#Simulatinng Measure-Z qubits

p00=basis(2,0) #initializaion of the data qubit
p01=basis(2,1)
pa=basis(2,0)
pb=basis(2,0)
pc=basis(2,0)
pd=basis(2,0)

psi_s=(p00+p01).unit()

t0=tensor(p00,pa,pb,pc,pd)
t1=tensor(p01,pa,pb,pc,pd)
t_s=tensor(psi_s,pa,pb,pc,pd)
print(t0)

#Measurement
def measure_comp_basis(psi):
    prob_0=abs(psi.overlap(t0))**2
    print(prob_0)
    prob_1=abs(psi.overlap(t1))**2
    print(prob_1)
    state=np.random.choice(2,p=[prob_0,prob_1])
    if(state==0):
        print("measurement: 0")
        psi=basis(2,0)
    elif(state==1):
        print("measurement: 1")
        psi=basis(2,1)
    return psi

f=t0.ptrace(0)
print(f)
print(f.data)

proj1=p00.proj()
print(proj1)
tmp1=p00+p01
print(tmp1)
tmp=proj1*tmp1
print(tmp)
