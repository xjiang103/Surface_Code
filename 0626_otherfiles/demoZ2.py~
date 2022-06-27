from IPython.display import Image
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
from qutip.qip.operations import cnot
from qutip.qip.circuit import QubitCircuit

#Simulatinng Measure-Z qubits

p00=basis(2,0) #initializaion of the data qubit
p01=basis(2,1)
psi_s=(p00+p01).unit()
t0=tensor(p00,p00,p00,p00,p00)
t1=tensor(p01,p00,p00,p00,p00)

t0_proj=tensor(p00.proj(),qeye(2),qeye(2),qeye(2),qeye(2))
t1_flip=tensor(sigmax(),qeye(2),qeye(2),qeye(2),qeye(2))
print(t0_proj)
dataqubit=[]
for i in range(4):
    rr=np.random.randint(0,2)
    print("random number:")
    print(rr)
    qtmp=(np.sqrt(rr)*p00+np.sqrt(1-rr)*p01).unit()
    dataqubit.append(qtmp)
    print(qtmp)
t_init=tensor(p00,dataqubit[0],dataqubit[1],dataqubit[2],dataqubit[3])
print(t_init)
print(t_init.ptrace(0))
print(t1.ptrace(0))
#Measurement
def measure_and_flip(psi):
    prob_0=abs((psi.ptrace(0)*p00).overlap(p00))**2
    prob_0_proj=prob_0
    prob_1=abs((psi.ptrace(0)*p01).overlap(p01))**2
    print("measure qubit:p0="+str(prob_0)+", p1="+str(prob_1))
    state=np.random.choice(2,p=[prob_0,prob_1])
    
    prob_0=(psi.ptrace(1)*p00).overlap(p00)
    prob_1=(psi.ptrace(1)*p01).overlap(p01)
    print("data qubit1:p0="+str(prob_0)+", p1="+str(prob_1))
    
    prob_0=(psi.ptrace(2)*p00).overlap(p00)
    prob_1=(psi.ptrace(2)*p01).overlap(p01)
    print("data qubit2:p0="+str(prob_0)+", p1="+str(prob_1))

    prob_0=(psi.ptrace(3)*p00).overlap(p00)
    prob_1=(psi.ptrace(3)*p01).overlap(p01)
    print("data qubit3:p0="+str(prob_0)+", p1="+str(prob_1))

    prob_0=(psi.ptrace(4)*p00).overlap(p00)
    prob_1=(psi.ptrace(4)*p01).overlap(p01)
    print("data qubit4:p0="+str(prob_0)+", p1="+str(prob_1))
    
    psitmp=psi
    if(prob_0_proj==0):
        psitmp=(t1_flip*(psi)).unit()
        return psitmp
    elif(state==0):
        print("measurement=0")
        psitmp=(t0_proj*(psi)).unit()
        return psitmp
    elif(state==1):
        print("measurement=1")
        psitmp=(t0_proj*(psi)).unit()
        return psitmp
#operation for stablization
def fourcnot(psi):
    q=QubitCircuit(5,reverse_states=False)
    q.add_gate("CNOT",controls=[1],targets=[0])
    q.add_gate("CNOT",controls=[2],targets=[0])
    q.add_gate("CNOT",controls=[3],targets=[0])
    q.add_gate("CNOT",controls=[4],targets=[0])
    qmat=gate_sequence_product(q.propagators())
    return qmat*psi

tstate=t_init
print(t_init.overlap(t0))
for count in range(10):
    print(str(count)+"th run")
    tstate=fourcnot(tstate)
    tstate=measure_and_flip(tstate)
