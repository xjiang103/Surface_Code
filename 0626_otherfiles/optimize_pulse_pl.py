from scipy.optimize import minimize
import subprocess
import numpy as np
import uuid
import os
from qutip import *
from qutip.qip.operations import snot,phasegate
import argparse

b=600
parser = argparse.ArgumentParser(description='Optimize Pulses')
parser.add_argument('-b','--b_couple', help='Rydberg Coupling Strength', required=False, type=float)
args = vars(parser.parse_args())
if args["b_couple"] is not None:
    b = args["b_couple"]

unique_file = str(uuid.uuid4())[0:8]
file_name = "dm_"+unique_file+".dat" #Allow us to run in parallel


def fun_sp(params,final_run=None):


    #Run QuaC
    try:
        output = subprocess.check_output(["./na_2_atom","-ts_rk_type","5bs","-ts_rtol","1e-8","-ts_atol","1e-8","-n_ens","-1",
                                          "-pulse_type","ARP","-file",file_name,
                                          "-b_term",str(b),
                                          "-delta",str(params[0]),
                                          "-pulse_length",str(params[1])])
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
    return 1-fid

def qutip_phase(params,dm):
    #define cz_arp
    cz_arp = Qobj([[1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,-1]],dims=[[2,2],[2,2]])

    #Get two hadamard state
    state = tensor(snot(),snot())*tensor(basis(2,0),basis(2,0))

    #Apply phase gates with parameters that we are optimizing
    state = tensor(qeye(2),phasegate(params[0]))*state
    state = tensor(phasegate(params[1]),qeye(2))*state

    #Now apply cz_arp
    state = cz_arp*state
    print("trace=",str(np.trace(dm)))
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
print("Final Params: ",str(res.x))
#Final Fidelity:  0.9997463238664505
#Final Fidelity:  0.9997626628800216
