from qutip import *
from qutip.qip.operations import cnot,snot

import numpy as np

#order of enumeration: 0 d 1 r
zero=basis(4,0)
one =basis(4,2)

#parameters

b_1r = 1.0/16.0
b_0r = 1.0/16.0
b_dr = 7.0/8.0

gamma_r = 1.0/(540.0)

b_field = 100000*(2*np.pi)

omega_amp=17.0*(2*np.pi)
delta_amp=23.0*(2*np.pi)
lpulse=0.54

#Hadamard Gate
def hada4():
    mat=np.array([[np.sqrt(2)/2,0,np.sqrt(2)/2,0],
                  [0,1,0,0],
                  [np.sqrt(2)/2,0,-1*np.sqrt(2)/2,0],
                  [0,0,0,1]])
    return Qobj(mat,dims=[[4],[4]])

#omega(t)
def H1_coeff(t,args):
    t1=lpulse/4.0
    t2=3*lpulse/4.0
    tau=0.175*lpulse
    a=np.exp(-1*np.power(t1,4)/np.power(tau,4))
    p=0.0
    if ((t>0)&(t<lpulse/2)):
        p=np.exp(-1*np.power(t-t1,4)/np.power(tau,4))/(1-a)
    elif ((t>lpulse/2)&(t<lpulse)):
        p=np.exp(-1*np.power(t-t2,4)/np.power(tau,4))/(1-a)
    else:
        p=0.0
    return omega_amp*p

#delta(t)
def H2_coeff(t,args):
    p=0.0
    if ((t>0)&(t<lpulse/2)):
        p=np.sin(2*np.pi*(t+3*lpulse/4)/lpulse)
    elif ((t>lpulse/2)&(t<lpulse)):
        p=np.sin(2*np.pi*(t+1*lpulse/4)/lpulse)
    else:
        p=0.0

    return delta_amp*p

#some operators for constructing the Hamiltonian and Lindblad term

mat_r1=np.array([[0,0,0,0],
                 [0,0,0,0],
                 [0,0,0,0],
                 [0,0,1,0]])
r_one=Qobj(mat_r1,dims=[[4],[4]])

mat_0r=np.array([[0,0,0,1],
                 [0,0,0,0],
                 [0,0,0,0],
                 [0,0,0,0]])
zero_r=Qobj(mat_0r,dims=[[4],[4]])

mat_1r=np.array([[0,0,0,0],
                 [0,0,0,0],
                 [0,0,0,1],
                 [0,0,0,0]])
one_r=Qobj(mat_1r,dims=[[4],[4]])

mat_dr=np.array([[0,0,0,0],
                 [0,0,0,1],
                 [0,0,0,0],
                 [0,0,0,0]])
d_r=Qobj(mat_dr,dims=[[4],[4]])

mat_rr=np.array([[0,0,0,0],
                 [0,0,0,0],
                 [0,0,0,0],
                 [0,0,0,1]]) #RR?
r_r=Qobj(mat_rr,dims=[[4],[4]])

#Hamiltonian and Lindblad term

hct=(1/2)*(r_one+one_r)

H0=b_field*tensor(r_r,r_r)
H1=tensor(hct,qeye(4))+tensor(qeye(4),hct)
H2=tensor(r_r,qeye(4))+tensor(qeye(4),r_r)

H=[H0,[H1,H1_coeff],[H2,H2_coeff]]

#Need to have lists of Linblads
c_ops = [np.sqrt(b_0r*gamma_r)*tensor(zero_r,qeye(4)),
         np.sqrt(b_1r*gamma_r)*tensor(one_r,qeye(4)),
         np.sqrt(b_dr*gamma_r)*tensor(d_r,qeye(4)),
         np.sqrt(b_0r*gamma_r)*tensor(qeye(4),zero_r),
         np.sqrt(b_1r*gamma_r)*tensor(qeye(4),one_r),
         np.sqrt(b_dr*gamma_r)*tensor(qeye(4),d_r)]

####
#Expectation values
#


#-------------------------------------------------------------------------------------------

data=tensor(one,one)
print("Initial State:")
print(data)


hada_1=tensor(hada4(),hada4())

data=hada_1*data
print("state after 1st Hadamard:")
print(data)

####
##apply cz gate
times = np.linspace(0.0, 0.54, 100000)


result = mesolve(H, data, times, c_ops, [])
####

print("State after CZ_ARP")
print(result.states[-1])

hada_2=tensor(qeye(4),hada4())
dm_hada2 = hada_2 * result.states[-1] * hada_2
print("state after 2nd Hadamard:")
print(dm_hada2)

print("two-level density matrix")
enum2=[0,2,8,10]
mat=np.zeros((4,4),dtype=complex)
for i in range(4):
    for j in range(4):
        mat[i][j]=dm_hada2[enum2[i]][0][enum2[j]]

twolvl_dm=Qobj(mat,dims=[[2,2],[2,2]])
print(twolvl_dm)


cz_arp = np.diag([1,-1,-1,-1])
psi0 = tensor(basis(2,1),basis(2,1))
psi0 = tensor(snot(),snot())*psi0
psi0 = cz_arp*psi0
psi0 = tensor(qeye(2),snot())*psi0
psi0 = ket2dm(Qobj(psi0))
psi0=Qobj(psi0,dims=[[2,2],[2,2]])

print(psi0)

print(fidelity(psi0,twolvl_dm))
