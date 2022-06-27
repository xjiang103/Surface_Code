from IPython.display import Image
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
from qutip.qip.operations import cnot,snot,swap
from qutip.qip.circuit import QubitCircuit
from qutip.qip.gates import gate_sequence_product
#Simulating Measure-Z qubits test

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

#Measurement and set the measurement(center) qubit back to ground state
def measure_and_flip(psi):
    #measure the measurement qubit
    print("measurement on the ancilla qubit:")
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
        psi_tmp = tensor(sigmax(),qeye(2),qeye(2),qeye(2),qeye(2))*psi_tmp #Reset ancilla to 0

#    psi_tmp=zero #initializaion for next round of operation: set the measurement qubit bact to ground state
    #psi_tmp=tensor(zero,dataqubit[0],dataqubit[1],dataqubit[2],dataqubit[3])

    return psi_tmp

#operation for stablization(4 CNOTs)
def user_gate1():
    diag=np.array([1,-1,-1,-1])
    mat=np.diag(diag)
    return Qobj(mat,dims=[[2,2],[2,2]])
cz_arp=user_gate1()
def fourczarp(psi):
    psi1=psi
    psi1=tensor(snot(),qeye(2),qeye(2),qeye(2),qeye(2))*psi1
    
    psi1=tensor(cz_arp,qeye(2),qeye(2),qeye(2))*psi1
    
    psi1=swap(5, targets=[1,2])*psi1
    psi1=tensor(cz_arp,qeye(2),qeye(2),qeye(2))*psi1
    psi1=swap(5, targets=[1,2])*psi1

    psi1=swap(5, targets=[1,3])*psi1
    psi1=tensor(cz_arp,qeye(2),qeye(2),qeye(2))*psi1
    psi1=swap(5, targets=[1,3])*psi1

    psi1=swap(5, targets=[1,4])*psi1
    psi1=tensor(cz_arp,qeye(2),qeye(2),qeye(2))*psi1
    psi1=swap(5, targets=[1,4])*psi1
    
    psi1=tensor(snot(),qeye(2),qeye(2),qeye(2),qeye(2))*psi1
    return psi1
        
def fourcnot(psi):
    q=QubitCircuit(5,reverse_states=False)
    q.add_gate("CNOT",controls=[1],targets=[0])
    q.add_gate("CNOT",controls=[2],targets=[0])
    q.add_gate("CNOT",controls=[3],targets=[0])
    q.add_gate("CNOT",controls=[4],targets=[0])
    qmat=gate_sequence_product(q.propagators())
    return qmat*psi

def state_print(psi):

    print("ancila qubit:")
    ptrace=psi.ptrace(0)
    prob_0=(zero_proj_dag*zero_proj*ptrace).tr()
    prob_1=(one_proj_dag*one_proj*ptrace).tr()
    print("prob_0=",str(prob_0),", prob_1=",str(prob_1))
    for i in range(4):
        print(make_ordinal(i+1)+" data qubit")
        ptrace=psi.ptrace(i+1)
        prob_0=(zero_proj_dag*zero_proj*ptrace).tr()
        prob_1=(one_proj_dag*one_proj*ptrace).tr()
        print("prob_0=",str(prob_0),", prob_1=",str(prob_1))
        state=np.random.choice(2,p=[prob_0,prob_1])
    print("------------------------------------------")
    return


    
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
print("initial state:")
state_print(tstate_init)
#Watch t0L stay in its state
print("prob_i=probability for the qubit to be in state |i>")
for count in range(10):
    tstate_prev = tstate
    print("--------------------------------------")
    print(str(count+1)+"th run")
    tstate=fourczarp(tstate)
    tstate=measure_and_flip(tstate)
    print("Fidelity with init: ",fidelity(tstate,tstate_init))
    print("Fidelity with previous: ",fidelity(tstate,tstate_prev))

print("-----------------------------------")
print("system is stablized at this final state:")
state_print(tstate)

print("=============================================================================")
############################################
# Encode data in 0L + 1L
############################################

#Get quiescent states

#initial state of the plaquette
t_init0L = tensor(zero,zero,zero,zero,zero)
t_init1L = tensor(zero,snot()*zero,snot()*zero,snot()*zero,snot()*zero)#Only data qubits in |+>

#Note that there can be multiple quiescent state, but once it collapses into one, it will state there
t1L = fourcnot(t_init1L)
t1L = measure_and_flip(t1L)

#operation on the plaquette

#Watch t1L stay in its state
print("Watch t1L stay in its state")
tstate=t1L
for count in range(5):
    print("--------------------------------------")
    print(str(count+1)+"th run")
    tstate=fourcnot(tstate)
    tstate=measure_and_flip(tstate)

#Watch t0L stay in its state
t0L = fourcnot(t_init0L)
t0L = measure_and_flip(t0L)

#Watch t0L stay in its state
print("Watch t1L stay in its state")
tstate=t0L
for count in range(6):
    print("--------------------------------------")
    print(str(count+1)+"th run")
    tstate=fourcnot(tstate)
    tstate=measure_and_flip(tstate)

