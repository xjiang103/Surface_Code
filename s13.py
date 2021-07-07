#simulation on the circuit in Fig 3b of https://arxiv.org/abs/1404.3747

from IPython.display import Image
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
from qutip.qip.operations import cnot,snot
from qutip.qip.circuit import QubitCircuit
from qutip.qip.gates import gate_sequence_product
#Simulating Measure-Z qubits test

#basis of qubit
zero=basis(2,0)
one=basis(2,1)
plus=snot()*basis(2,0)
minus=snot()*basis(2,1)


#density matrix and projection operator of the basis states
zero_mat=ket2dm(zero)
one_mat=ket2dm(one)
plus_mat=ket2dm(plus)
minus_mat=ket2dm(minus)

zero_proj=zero.proj()
zero_proj_dag=zero_proj.dag()
one_proj=one.proj()
one_proj_dag=one_proj.dag()

plus_proj=plus.proj()
plus_proj_dag=plus_proj.dag()
minus_proj=minus.proj()
minus_proj_dag=minus_proj.dag()

#Measurement and set the ancila(center) qubit back to ground state
def ancilla_measure_reset(psi):
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
        #Collapse the ancilla to 0
        project_ancilla = tensor(zero_proj,qeye(2),qeye(2),qeye(2),qeye(2))
        psi_tmp = project_ancilla*psi
        psi_tmp = psi_tmp.unit() # Renormalize!
    elif(state==1):
        print("measurment=1")
        #collapse the ancilla to 1, reset to 0
        project_ancilla = tensor(one_proj,qeye(2),qeye(2),qeye(2),qeye(2))
        psi_tmp = project_ancilla*psi
        psi_tmp = psi_tmp.unit() # Renormalize!
        psi_tmp = tensor(sigmax(),qeye(2),qeye(2),qeye(2),qeye(2))*psi_tmp

    return psi_tmp

def ancilla_hadamard(psi):
    return tensor(snot(),qeye(2),qeye(2),qeye(2),qeye(2))*psi

#operation for stablization(4 CNOTs)
def xcnot(psi):
    q=QubitCircuit(5,reverse_states=False)
    q.add_gate("CNOT",controls=[0],targets=[3])
    q.add_gate("CNOT",controls=[0],targets=[4])
    q.add_gate("CNOT",controls=[0],targets=[1])
    q.add_gate("CNOT",controls=[0],targets=[2])
    qmat=gate_sequence_product(q.propagators())
    return qmat*psi

def zcnot(psi):
    q=QubitCircuit(5,reverse_states=False)
    q.add_gate("CNOT",controls=[1],targets=[0])
    q.add_gate("CNOT",controls=[3],targets=[0])
    qmat=gate_sequence_product(q.propagators())
    return qmat*psi

def make_ordinal(n):
    n = int(n)
    suffix = ['th', 'st', 'nd', 'rd', 'th'][min(n % 10, 4)]
    if 11 <= (n % 100) <= 13:
        suffix = 'th'
    return str(n) + suffix
######################################################
#Watch unstabilized random state become stabilized
#######################################################
data = rand_ket(16) #Random entangled state of the data
data.dims = [[2,2,2,2],[1,1,1,1]] #Set the dimensions so it looks like 4 qubits
tstate_init = tensor(zero,data)
tstate = tensor(zero,data)
print("prob_i=probability for the qubit to be in state |i>")
#Watch t0L stay in its state
for count in range(5):
    print("--------------------------------------")
    print(make_ordinal(count+1)+" run")
    
    tstate_prev1 = tstate    
    tstate=ancilla_hadamard(tstate)
    tstate=xcnot(tstate)
    tstate=ancilla_hadamard(tstate)
    print("1st measurement")
    tstate=ancilla_measure_reset(tstate)
    print("Fidelity with init: ",fidelity(tstate,tstate_init))
    print("Fidelity with previous: ",fidelity(tstate,tstate_prev1))

    tstate_prev2 = tstate   
    tstate=zcnot(tstate)
    print("2nd measurement:")
    tstate=ancilla_measure_reset(tstate)
    print("Fidelity with init: ",fidelity(tstate,tstate_init))
    print("Fidelity with previous: ",fidelity(tstate,tstate_prev2))
    

