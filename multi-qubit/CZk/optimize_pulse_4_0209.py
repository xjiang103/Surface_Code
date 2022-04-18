from scipy.optimize import minimize
import subprocess
import numpy as np
import uuid
import os
from qutip import *
from qutip.qip.operations import snot,phasegate
import argparse


parser = argparse.ArgumentParser(description='Optimize Pulses')
parser.add_argument('-b','--b_couple', help='Rydberg Coupling Strength', required=True, type=float)
args = vars(parser.parse_args())
b = args["b_couple"]
ddfac=0.2
typestr="cz3"


unique_file = str(uuid.uuid4())[0:8]
file_name = "dm_"+unique_file+".dat" #Allow us to run in parallel
params = [-0.5,0.2]

f=open("op_4_par_0119.txt","a")
f.write('\n')
f.write(typestr+' ')
def fun_sp(params):


    #Run QuaC
    try:
        output = subprocess.check_output(["./na_4_par_3lvl2_0209","-ts_rk_type","5bs","-ts_rtol","1e-8","-ts_atol","1e-8","-n_ens","-1",
                                          "-pulse_type","SP","-file",file_name,
                                          #"-bitstr","1111"
                                          "-b_term",str(params[2]),
                                          "-delta",str(params[0]),
                                          "-deltat",str(params[1]),
                                          "-dd_fac",str(ddfac)])
    except:
        pass

    #Read in the QuaC DM
    dm = Qobj(np.loadtxt(file_name).view(complex),dims=[[2,2,2,2],[2,2,2,2]])
    #Remove file
    os.remove(file_name)

    #QUTIP to get perfect circuit
    res = qutip_phase([0.10896237,-1.64590753,-1.6458525,-1.64582669],dm)

    fid = 1-res
    print(fid)
##    if(final_run):
##    print("Phase: ",res.x)
##    f.write(str(res.x[0])+' ')
##    f.write(str(res.x[1])+' ')
##    f.write(str(res.x[2])+' ')
##    f.write(str(res.x[3])+' ')
    return 1-fid

def qutip_phase(params,dm):
    #define cz_arp and czz arp
    cccz_arp = Qobj(np.diag([1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])
                   ,dims=[[2,2,2,2],[2,2,2,2]])
    czzz_arp = Qobj(np.diag([1,-1,-1,1,-1,1,1,-1,-1,-1,-1,1,-1,1,1,-1])
                   ,dims=[[2,2,2,2],[2,2,2,2]])
    
    #Get two hadamard state
    state = tensor(basis(2,0),basis(2,1),basis(2,0),basis(2,0))

    #Apply phase gates with parameters that we are optimizing
    state = tensor(phasegate(params[0]),qeye(2),qeye(2),qeye(2))*state
    state = tensor(qeye(2),phasegate(params[1]),qeye(2),qeye(2))*state
    state = tensor(qeye(2),qeye(2),phasegate(params[2]),qeye(2))*state
    state = tensor(qeye(2),qeye(2),qeye(2),phasegate(params[3]))*state

    #Now apply cz_arp
    state = czzz_arp*state

    #Get fidelity wrt quac dm
    fid = fidelity(dm,state)

    return 1-fid

def fun_arp(delta):
    #NOT COMPLETED!
    return 1-fid


#print("Optimizing SP for b = ",str(b))
f.write(str(ddfac)+' ')
#f.write(str(b)+' ')
default_sp_params = [-0.5,0.2]
default_sp_params = [0.49790076,0.16799145,10.75]
fidf=1-fun_sp(default_sp_params)

print("fidelity="+str(fidf))
#get the optimal phases
##fun_sp(res.x,True)
##print("Final Fidelity: ",str(1-res.fun))
##f.write(str(1-res.fun)+' ')
##print("Final Params: ",str(res.x))
##f.write(str(res.x[0])+' ')
##f.write(str(res.x[1])+' ')
##f.write(str(res.x[2])+' ')
##f.write('\n')
###Final Fidelity:  0.9997463238664505
f.close()
