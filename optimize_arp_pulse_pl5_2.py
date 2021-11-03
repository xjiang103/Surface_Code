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
ddfac=1
typestr="c4z"
print("b="+str(b))

unique_file = str(uuid.uuid4())[0:8]
file_name = "dm_"+unique_file+".dat" #Allow us to run in parallel
params = [-0.5,0.2,10]


f=open("op_5_arp_par.txt","a")
f.write('\n')
f.write(typestr+' ')
output = subprocess.check_output(["./na_5_par_3lvl3","-ts_rk_type","5bs","-ts_rtol","1e-8","-ts_atol","1e-8","-n_ens","-1",
                                    "-pulse_type","ARP","-file",file_name,
                                    "-b_term",str(b),
                                    "-delta",str(params[0]),
                                    "-pulse_length",str(params[1]),
                                    "-dd_fac",str(ddfac)])

    #Read in the QuaC DM
dm = Qobj(np.loadtxt(file_name).view(complex),dims=[[2,2,2,2,2],[2,2,2,2,2]])
    #Remove file


    #QUTIP to get perfect circuit
c4z_arp = Qobj(np.diag([1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])
                ,dims=[[2,2,2,2,2],[2,2,2,2,2]])
cz4_arp = Qobj(np.diag([1,-1,-1,1,-1,1,1,-1,-1,1,1,-1,1,-1,-1,1,-1,-1,-1,1,-1,1,1,-1,-1,1,1,-1,1,-1,-1,1])
                ,dims=[[2,2,2,2,2],[2,2,2,2,2]])
    
    #Get two hadamard state
state = tensor(snot(),snot(),snot(),snot(),snot())*tensor(basis(2,1),basis(2,1),basis(2,1),basis(2,1),basis(2,1))
params=[0,0,0,0,0]
    #Apply phase gates with parameters that we are optimizing
state = tensor(phasegate(params[0]),qeye(2),qeye(2),qeye(2),qeye(2))*state
state = tensor(qeye(2),phasegate(params[1]),qeye(2),qeye(2),qeye(2))*state
state = tensor(qeye(2),qeye(2),phasegate(params[2]),qeye(2),qeye(2))*state
state = tensor(qeye(2),qeye(2),qeye(2),phasegate(params[3]),qeye(2))*state
state = tensor(qeye(2),qeye(2),qeye(2),qeye(2),phasegate(params[4]))*state

    #Now apply cz_arp
state = c4z_arp*state
print("trace=",str(np.trace(dm)))
    #Get fidelity wrt quac dm
fid = fidelity(dm,state)



print("Optimizing ARP for b = ",str(b))
print("Optimizing Delta, T, and phases")


print("Final Fidelity: ",str(fid))
f.write(str(fid))

f.write('\n')
#Final Fidelity:  0.9997463238664505
f.close()
