from qutip import *
from qutip.qip.operations import cnot,snot,phasegate
import numpy as np
from scipy.optimize import minimize


f=open("quac_mat_2.txt",'r')

quac_mat=np.zeros((81,81),dtype=complex)
for i in range(19):
    r=f.readline()
    j,k,real,imag=r.split()
    quac_mat[int(j)][int(k)]=int(real)+1j*int(imag)
quac_m=Qobj(quac_mat,dims=[[3,3,3,3],[3,3,3,3]])
#order of enumeration: 0 1 r
zero=basis(3,0)
one =basis(3,1)
r   =basis(3,2)

#parameters

b_1r = 1.0/16.0
b_0r = 1.0/16.0
b_dr = 7.0/8.0

gamma_r = 0/(540.0)

b_field = 600*(2*np.pi)

omega_amp=17.0*(2*np.pi)

delta_200=23.85969*(2*np.pi)
delta_300=23.94195557*(2*np.pi)
delta_400=24.01839294*(2*np.pi)
delta_500=23.86193848*(2*np.pi)
delta_600=23.86193848*(2*np.pi)
delta_amp=delta_200

lpulse=0.54

params_200=[0.57267,0.57264,0.57279,0.57249]
params_300=[0.38050708,0.38048573,0.3802977,0.38054245] 
params_400=[0.2848848,0.28477299,0.28486211,0.28502494]
params_500=[0.22773865,0.22759277,0.22770762,0.22759728] 
params=params_200
#Hadamard Gate
def hada3():
    mat=np.array([[np.sqrt(2)/2,np.sqrt(2)/2,0],
                  [np.sqrt(2)/2,-1*np.sqrt(2)/2,0],
                  [0,0,1]])
    return Qobj(mat,dims=[[3],[3]])

#omega(t)
def Omega_coeff(t,args):
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
def Delta_coeff(t,args):
    p=0.0
    if ((t>0)&(t<lpulse/2)):
        p=np.sin(2*np.pi*(t+3*lpulse/4)/lpulse)
    elif ((t>lpulse/2)&(t<lpulse)):
        p=np.sin(2*np.pi*(t+1*lpulse/4)/lpulse)
    else:
        p=0.0

    return delta_amp*p

def qutip_phase(params,qvec):
    #define cz_arp and czz arp
    cccz_arp = Qobj(np.diag([1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])
                   ,dims=[[2,2,2,2],[2,2,2,2]])
    czzz_arp = Qobj(np.diag([1,-1,-1,1,-1,1,1,-1,-1,-1,-1,1,-1,1,1,-1])
                   ,dims=[[2,2,2,2],[2,2,2,2]])
    
    #Get two hadamard state
    state = tensor(snot(),snot(),snot(),snot())*tensor(basis(2,1),basis(2,1),basis(2,1),basis(2,1))

    #Apply phase gates with parameters that we are optimizing
    state = tensor(phasegate(params[0]),qeye(2),qeye(2),qeye(2))*state
    state = tensor(qeye(2),phasegate(params[1]),qeye(2),qeye(2))*state
    state = tensor(qeye(2),qeye(2),phasegate(params[2]),qeye(2))*state
    state = tensor(qeye(2),qeye(2),qeye(2),phasegate(params[3]))*state

    #Now apply cz_arp
    state = cccz_arp*state
    #print("trace=",str(np.trace(dm)))

    #Get fidelity wrt quac dm
    fid = fidelity(qvec,state)

    return 1-fid
    
#some operators for constructing the Hamiltonian and Lindblad term
r_one=r*one.dag()

zero_r=zero*r.dag()

one_r=one*r.dag()

r_r=r*r.dag()

#Hamiltonian and Lindblad term

hct=(1/2)*(r_one+one_r)

#assuming the same coupling for every pair of atoms

H0=b_field*(tensor(r_r,r_r,qeye(3),qeye(3))+
            tensor(r_r,qeye(3),r_r,qeye(3))+
            tensor(r_r,qeye(3),qeye(3),r_r)+
            tensor(qeye(3),r_r,r_r,qeye(3))+
            tensor(qeye(3),r_r,qeye(3),r_r)+
            tensor(qeye(3),qeye(3),r_r,r_r))
##H0=b_field*(tensor(r_r,r_r,qeye(3),qeye(3))+
##            tensor(r_r,qeye(3),r_r,qeye(3))+
##            tensor(r_r,qeye(3),qeye(3),r_r))

H_leak=(-1/1)*(1j)*b_dr*gamma_r*(tensor(r_r,qeye(3),qeye(3),qeye(3))+
                               tensor(qeye(3),r_r,qeye(3),qeye(3))+
                               tensor(qeye(3),qeye(3),r_r,qeye(3))+
                               tensor(qeye(3),qeye(3),qeye(3),r_r))

H1=(tensor(hct,qeye(3),qeye(3),qeye(3))+
    tensor(qeye(3),hct,qeye(3),qeye(3))+
    tensor(qeye(3),qeye(3),hct,qeye(3))+
    tensor(qeye(3),qeye(3),qeye(3),hct))

H2=(tensor(r_r,qeye(3),qeye(3),qeye(3))+
    tensor(qeye(3),r_r,qeye(3),qeye(3))+
    tensor(qeye(3),qeye(3),r_r,qeye(3))+
    tensor(qeye(3),qeye(3),qeye(3),r_r))

H=[H0+H_leak,[H1,Omega_coeff],[H2,Delta_coeff]]

print("H=:\n")
print(H1.data)

#---------------------------------------------------------------------

data=tensor(one,one,one,one)

hada_1=tensor(hada3(),hada3(),hada3(),hada3())
data=hada_1*data

####
#apply C_3 Z gate
times=np.linspace(0.0,0.54,10000)

options = Options(normalize_output=False)
result = sesolve(H,data,times,options=options)
##result=sesolve(H,data,times)

res=result.states[-1]

print("two-level density matrix")
print(res)
print(res.norm())

enum2=[0,1,3,4,9,10,12,13,27,28,30,31,36,37,39,40]
qvecmat=np.zeros((16),dtype=complex)
for i in range(16):
    qvecmat[i]=res[enum2[i]]
qvec=Qobj(qvecmat,dims=[[2,2,2,2],[1,1,1,1]])
print(qvec)
print(qvec.norm())

res = minimize(qutip_phase,[0,0,0,0],method="COBYLA",args=(qvec))
fid = 1-res.fun
print("Fidelity=",fid)
print("Phase: ",res.x)
##enum2=[
##dmmat=np.zeros((16,16),dtype=complex)
##for i in range



