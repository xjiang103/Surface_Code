from IPython.display import Image
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
from qutip.qip.operations import cnot
from qutip.qip.circuit import QubitCircuit
from qutip.qip.gates import gate_sequence_product
#Simulating Measure-Z qubits

#basis of qubit
zero=basis(2,0) 
one=basis(2,1)

#density matrix and projection operator of the basis states
zero_mat=ket2dm(zero)
one_mat=ket2dm(one)

zero_proj=zero.proj()
zero_proj_dag=zero_proj.dag()
one_proj=one.proj()
one_proj_dag=one_proj.dag()

#Initialization: the data qubit is initially in a arbitrary entangled state
dataqubit=[]
print("Initialization:")
for i in range(4):
    rr=np.random.random()
    datatmp=(np.sqrt(rr)*zero+np.sqrt(1-rr)*one).unit()
    dataqubit.append(datatmp)
    print(str(i)+"th data qubit initially at:"+str(np.sqrt(rr))+"|0> + "+str(np.sqrt(1-rr))+"|1>")

#initial state of the plaquette
t_init=tensor(zero,dataqubit[0],dataqubit[1],dataqubit[2],dataqubit[3])

#Measurement and set the measurement(center) qubit back to ground state
def measure_and_flip(psi):
    #measure the measurement qubit
    ptrace=psi.ptrace(0)

    #calculate the probability of obtaining each result
    prob_0=(zero_proj_dag*zero_proj*ptrace).tr()
    prob_1=(one_proj_dag*one_proj*ptrace).tr()
    print("prob_0=",str(prob_0),", prob_1=",str(prob_1))

    #measurement result
    state=np.random.choice(2,p=[prob_0,prob_1])
    
    if(state==0):
        print("measurement=0")
    elif(state==1):
        print("measuremnt=1")
        
    psi_tmp=zero #initializaion for next round of operation: set the measurement qubit bact to ground state

    for i in range(4):
        #measure the data qubit
        ptrace=psi.ptrace(i+1)

        #calculate the probability of obtaining each result
        prob_0=(zero_proj_dag*zero_proj*ptrace).tr()
        prob_1=(one_proj_dag*one_proj*ptrace).tr()

        #measurement result
        state=np.random.choice(2,p=[prob_0,prob_1])

        #this qubit collapses to the state that corresponds to the measurement result 
        if(state==0):
            psi_tmp=tensor(psi_tmp,zero)
            print(str(i+1)+"th data qubit in state |"+str(state)+">")
        elif(state==1):
            psi_tmp=tensor(psi_tmp,one)
            print(str(i+1)+"th data qubit in state |"+str(state)+">")
    return psi_tmp

#operation for stablization(4 CNOTs)
def fourcnot(psi):
    q=QubitCircuit(5,reverse_states=False)
    q.add_gate("CNOT",controls=[1],targets=[0])
    q.add_gate("CNOT",controls=[2],targets=[0])
    q.add_gate("CNOT",controls=[3],targets=[0])
    q.add_gate("CNOT",controls=[4],targets=[0])
    qmat=gate_sequence_product(q.propagators())
    return qmat*psi

#operation on the plaquette
tstate=t_init
for count in range(6):
    print("--------------------------------------")
    print(str(count+1)+"th run")
    tstate=fourcnot(tstate)
    tstate=measure_and_flip(tstate)
