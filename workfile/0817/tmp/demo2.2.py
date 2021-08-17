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
print(t_s)


t_s.eliminate_states([0,1,2,3])
print(t_s)

td=tensor(p01,t_s)
print(t_s)
