from scipy.optimize import minimize
import subprocess
import numpy as np
import uuid
import os
from qutip import *
from qutip.qip.operations import snot,phasegate
import argparse

#enter "xxx" for the equal superposition state. 
parser = argparse.ArgumentParser(description='Optimize Pulses')
parser.add_argument('-b','--b_couple', help='Rydberg Coupling Strength', required=True, type=float)
parser.add_argument('-ist','--initial_state', help='Variation in Rydberg Coupling Strength', required=True, type=int)
args = vars(parser.parse_args())
b = args["b_couple"]
ist=args["initial_state"]
print("initial state is "+str(ist))
ddfac=1
typestr="ccz"

unique_file = str(uuid.uuid4())[0:8]
file_name = "dm_"+unique_file+".dat" #Allow us to run in parallel

f=open("arp_2_F.txt","a")
f.write('\n')
f.write("initial state="+str(ist)+'\n')
f.write(typestr+' ')

def ia(init_state):
    state_arr=[]
    for k in range(len(init_state)):
        if (init_state[k]=='0'):
            state_arr.append(basis(2,0))
        elif (init_state[k]=='1'):
            state_arr.append(basis(2,1))
    state=tensor(state_arr[0],state_arr[1])
    return state
init_arr=[]
str_arr=[]

init_arr.append(ia("00"))
init_arr.append(ia("01"))
init_arr.append(ia("10"))
init_arr.append(ia("11"))
ave_state=(1.0/np.sqrt(4))*(init_arr[0]+init_arr[1]+init_arr[2]+init_arr[3])
#ave_state=(1.0/np.sqrt(2))*(init_arr[0]+init_arr[1])
init_arr.append(ave_state)

str_arr=["00","01","10","11","xx"]
#str_arr=["00","11","xx"]
init_state=init_arr[ist]
print("initial state is "+str_arr[ist])
def fun_sp(params,final_run=None):


    #Run QuaC
    try:
        output = subprocess.check_output(["./na_1115","-ts_rk_type","5bs","-ts_rtol","1e-8","-ts_atol","1e-8","-n_ens","-1",
                                          "-pulse_type","ARP","-file",file_name,
                                          "-bitstr",str_arr[ist],
                                          "-b_term",str(b),
                                          "-delta",str(params[0]),
                                          "-pulse_length",str(params[1]),
                                          "-dd_fac",str(ddfac)])
    except:
        pass

    #Read in the QuaC DM
    dm = Qobj(np.loadtxt(file_name).view(complex),dims=[[2,2],[2,2]])
    #Remove file
    os.remove(file_name)

    #QUTIP to get perfect circuit
    res = minimize(qutip_phase,[0,0],method="COBYLA",args=(dm))

    fid = 1-res.fun
    print(fid)
    if(final_run):
        print("Phase: ",res.x)
        f.write(str(res.x[0])+' ')
        f.write(str(res.x[1])+' ')
    return 1-fid

def qutip_phase(params,dm):
    #define cz_arp
    cz_arp = Qobj([[1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,-1]],dims=[[2,2],[2,2]])

    #Get two hadamard state
    #state = tensor(snot(),snot())*tensor(basis(2,0),basis(2,0))
    state=init_state

    #Apply phase gates with parameters that we are optimizing
    state = tensor(qeye(2),phasegate(params[0]))*state
    state = tensor(phasegate(params[1]),qeye(2))*state

    #Now apply cz_arp
    state = cz_arp*state
    #print("trace=",str(np.trace(dm)))
    #Get fidelity wrt quac dm
    fid = fidelity(dm,state)

    return 1-fid

def fun_arp(delta):
    #NOT COMPLETED!
    return 1-fid



print("Optimizing SP for b = ",str(b))

default_sp_params = [23,0.54]
res = minimize(fun_sp,default_sp_params,method="nelder-mead")

#get the optimal phases
fun_sp(res.x,True)
print("Final Fidelity: ",str(1-res.fun))
f.write(str(1-res.fun)+' ')
print("Final Params: ",str(res.x))
f.write(str(res.x[0])+' ')
f.write(str(res.x[1])+' ')
f.write('\n')
#Final Fidelity:  0.9997463238664505
#Final Fidelity:  0.9997626628800216
