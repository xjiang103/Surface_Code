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
parser.add_argument('-deltab','--deltab_couple', help='Variation in Rydberg Coupling Strength', required=True, type=float)
args = vars(parser.parse_args())
b = args["b_couple"]
deltab=args["deltab_couple"]

ddfac=1
typestr="ccz"


unique_file = str(uuid.uuid4())[0:8]
file_name = "dm_"+unique_file+".dat" #Allow us to run in parallel
params = [-0.5,0.2]


f=open("op_3_arp_par_0811.txt","a")
f.write('\n')
f.write(typestr+' ')
#params=[23.682658958435056,0.5521307086944582]
params=[23.754459381103526,0.5530808029174807]

subprocess.run(["./na_3_par_3lvl2_0809","-ts_rk_type","5bs","-ts_rtol","1e-8","-ts_atol","1e-8","-n_ens","-1",
                                          "-pulse_type","ARP","-file",file_name,
                                          "-b_term",str(b),
                                          "-delta_b_term",str(deltab),
                                          "-delta",str(params[0]),
                                          "-pulse_length",str(params[1]),
                                          "-dd_fac",str(ddfac)])

dm = Qobj(np.loadtxt(file_name).view(complex),dims=[[2,2,2],[2,2,2]])

ccz_arp = Qobj([[1,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0],[0,0,0,-1,0,0,0,0],[0,0,0,0,-1,0,0,0],[0,0,0,0,0,-1,0,0],[0,0,0,0,0,0,-1,0],[0,0,0,0,0,0,0,-1]],dims=[[2,2,2],[2,2,2]])
czz_arp = Qobj([[1,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1]],dims=[[2,2,2],[2,2,2]])
    
    #Get two hadamard state
state = tensor(snot(),snot(),snot())*tensor(basis(2,1),basis(2,1),basis(2,1))

#params=[0.49044988054774474,0.4903519642692983,0.4903430001222758]
params=[0.24433141136863257,0.2441608870373369,0.24423646035956137]

    #Apply phase gates with parameters that we are optimizing
state = tensor(phasegate(params[0]),qeye(2),qeye(2))*state
state = tensor(qeye(2),phasegate(params[1]),qeye(2))*state
state = tensor(qeye(2),qeye(2),phasegate(params[2]))*state

    #Now apply cz_arp
state = ccz_arp*state
print("trace=",str(np.trace(dm)))
    #Get fidelity wrt quac dm
fid = fidelity(dm,state)
print("trace=",str(np.trace(dm)))

print("Optimizing ARP for b = ",str(b))
print("Optimizing Delta, T, and phases")
f.write(str(ddfac)+' ')
#f.write(str(b)+' ')
f.write("Delta_T_phases for b="+str(b)+' ')


print("Final Fidelity: ",fid)
f.write(str(fid)+' ')

#Final Fidelity:  0.9997463238664505
f.close()
