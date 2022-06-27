from qutip import *
from qutip.qip.operations import cnot,snot,phasegate
import numpy as np
from scipy.optimize import minimize
from tabulate import tabulate
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 16})


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

b_field = 200*(2*np.pi)

omega_amp=17.0*(2*np.pi)

delta_200=23.85969238*(2*np.pi)
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
        p=(np.exp(-1*np.power(t-t1,4)/np.power(tau,4))-a)/(1-a)
    elif ((t>lpulse/2)&(t<lpulse)):
        p=(np.exp(-1*np.power(t-t2,4)/np.power(tau,4))-a)/(1-a)
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

#H=[H0+H_leak,[H1,Omega_coeff],[H2,Delta_coeff]]
H=[H0,[H1,Omega_coeff],[H2,Delta_coeff]]


#---------------------------------------------------------------------

data=tensor(one,one,one,one)

hada_1=tensor(hada3(),hada3(),hada3(),hada3())
data=hada_1*data

####
#apply C_3 Z gate

times=np.linspace(0.0,0.54,10000)
options = Options(normalize_output=False,atol=1e-15,rtol=1e-15,nsteps=10000,tidy=False)
result = mesolve(H,data,times,options=options)
##result=sesolve(H,data,times)

#res_mat=result.states

state1111=[]
state111r=[]
state11r1=[]
state1r11=[]
stater111=[]
residual=[]
enum2=[0,1,3,4,9,10,12,13,27,28,30,31,36,37,39,40]


print(result.states[0][13][0][0])
for i in range(10000):
    state1111.append(np.abs(result.states[i][40][0][0])**2)
    state111r.append(np.abs(result.states[i][41][0][0])**2)
    state11r1.append(np.abs(result.states[i][43][0][0])**2)
    state1r11.append(np.abs(result.states[i][49][0][0])**2)
    stater111.append(np.abs(result.states[i][67][0][0])**2)
    residual.append(1-state1111[i]-state111r[i]-state11r1[i]-state1r11[i]-stater111[i])


plt.plot(times+0.00,state1111,label=r"qutip $|1111\rangle\langle 1111|$")
plt.plot(times+0.00,state111r,label=r"qutip $|111r\rangle\langle 111r|$")
#plt.plot(times+0.00,state11r1,label=r"$qutip |11r1\rangle\langle 11r1|$")
#plt.plot(times+0.00,state1r11,label=r"$qutip |1r11\rangle\langle 1r11|$")
#plt.plot(times+0.00,stater111,label=r"$qutip |r111\rangle\langle r111|$")
plt.plot(times+0.00,residual,label="qutip other states")

quacdata = np.loadtxt("na_4_par_3lvl2_0714.txt")

plt.plot(quacdata[:,1],quacdata[:,2],lw=2,label=r"quac $|1111\rangle\langle 1111|$")
plt.plot(quacdata[:,1],quacdata[:,3],lw=2,label=r"quac $|111r\rangle\langle 111r|$")
#plt.plot(quacdata[:,1],quacdata[:,4],lw=2,label=r"quac $|11r1\rangle\langle 11r1|$")
#plt.plot(quacdata[:,1],quacdata[:,5],lw=2,label=r"quac $|1r11\rangle\langle 1r11|$")
#plt.plot(quacdata[:,1],quacdata[:,6],lw=2,label=r"quac $|r111\rangle\langle r111|$")
plt.plot(quacdata[:,1],1-quacdata[:,6]-quacdata[:,5]-quacdata[:,4]-quacdata[:,3]-quacdata[:,2],lw=2,label=r"quac other states")


##residual=[]
##enum2=[0,1,3,4,9,10,12,13,27,28,30,31,36,37,39,40]
##
##for i in range(10000):
##    data_tmp=0
##    for j in range(16):
##        data_tmp=data_tmp+np.abs(result.states[i][enum2[j]][0][0])**2
##    residual.append(1-data_tmp)
##plt.plot(times,residual,label="None-computational states")


plt.title("B= 200MHz")
plt.xlabel("Time")
plt.xlim(0,0.6)
plt.yscale("log")
plt.ylim(1e-7,1)
plt.legend(bbox_to_anchor=(0.54, 0.4), loc='center left', ncol=1)

plt.show()

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



