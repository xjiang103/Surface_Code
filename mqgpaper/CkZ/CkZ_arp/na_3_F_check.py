from scipy.optimize import minimize
import subprocess
import numpy as np
import uuid
import os
from qutip import *
from qutip.qip.operations import snot,phasegate
import argparse
from joblib import Parallel, delayed


parser = argparse.ArgumentParser(description='Optimize Pulses')
parser.add_argument('-b','--b_couple', help='Rydberg Coupling Strength', required=True, type=float)
#parser.add_argument('-ist','--initial_state', help='Variation in Rydberg Coupling Strength', required=True, type=int)
args = vars(parser.parse_args())
b = args["b_couple"]
#ist=args["initial_state"]
#print("initial state is "+str(ist))
ddfac=1
typestr="ccz"


unique_file = str(uuid.uuid4())[0:8]
file_name = "dm_"+unique_file+".dat" #Allow us to run in parallel
params = [23]

def ia(init_state):
    state_arr=[]
    for k in range(len(init_state)):
        if (init_state[k]=='0'):
            state_arr.append(basis(2,0))
        elif (init_state[k]=='1'):
            state_arr.append(basis(2,1))
    state=tensor(state_arr[0],state_arr[1],state_arr[2])
    return state
init_arr=[]
str_arr=[]

#init_arr is an arry of initial states
#str_arr is an arry of strings, to be used in calling the QuaC program

#enumerate the basis states
for j1 in range(2):
    for j2 in range(2):
        for j3 in range(2):
            initstr=""
            if (j1==0):
                initstr=initstr+'0'
            elif (j1==1):
                initstr=initstr+'1'
            if (j2==0):
                initstr=initstr+'0'
            elif (j2==1):
                initstr=initstr+'1'
            if (j3==0):
                initstr=initstr+'0'
            elif (j3==1):
                initstr=initstr+'1'
            str_arr.append(initstr)
            init_arr.append(ia(initstr))
ave_state=0

def qutip_phase(params,i,dm):
    #define cz_arp and czz arp
    ccz_arp = Qobj([[1,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0],[0,0,0,-1,0,0,0,0],[0,0,0,0,-1,0,0,0],[0,0,0,0,0,-1,0,0],[0,0,0,0,0,0,-1,0],[0,0,0,0,0,0,0,-1]],dims=[[2,2,2],[2,2,2]])
    czz_arp = Qobj([[1,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1]],dims=[[2,2,2],[2,2,2]])
    fid=0
    state=init_arr[i]
    #Apply phase gates with parameters that we are optimizing
    state = tensor(phasegate(params[0]),qeye(2),qeye(2))*state
    state = tensor(qeye(2),phasegate(params[1]),qeye(2))*state
    state = tensor(qeye(2),qeye(2),phasegate(params[2]))*state

    #Now apply cz_arp
    state = ccz_arp*state

    #Get fidelity wrt quac dm
    fid_tmp = fidelity(dm,state)
    return fid_tmp

#add the additional superposition state
for j in range(8):
    ave_state+=(1/np.sqrt(8))*init_arr[j]
init_arr.append(ave_state)
str_arr.append("xxx")
params = [23.824314452707778, 0.5353643417805434 ]
for i in range(9):
    unique_file = str(uuid.uuid4())[0:8]
    file_name = "dm_"+unique_file+".dat" #Allow us to run in parallel
    #Run QuaC
    try:
        output = subprocess.check_output(["./na_3_F_par","-ts_rk_type","5bs","-ts_rtol","1e-8","-ts_atol","1e-8","-n_ens","-1",
                                          "-pulse_type","ARP","-file",file_name,
                                          "-bitstr",str_arr[i],
                                          "-b_term",str(b),
                                          "-delta",str(params[0]),
                                          "-pulse_length",str(params[1]),
                                          "-dd_fac",str(ddfac)])
    except:
        pass

    #Read in the QuaC DM

    dm = Qobj(np.loadtxt(file_name).view(complex),dims=[[2,2,2],[2,2,2]])
    #Remove file
    os.remove(file_name)
    print("fid ",i,": ",qutip_phase([0.49034015738225945, 0.4905471860442265, 0.49051937108012406 ],i,dm))

