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
typestr="ccz"

fnarr=[]
for i in range(8):  
    unique_file = str(uuid.uuid4())[0:8]
    file_name = "dm_"+unique_file+".dat" #Allow us to run in parallel
    fnarr.append(file_name)
init_state=["000","001","010","011","100","101","110","111"]

f=open("op_3_arp_par_F2.txt","a")
f.write('\n')
f.write(typestr+' ')
def fun_sp(params,final_run=None):
    fidarr=[]
    earr=[]
    for i in range(8):
        print(i)
        print(init_state[i])
        try:
            output = subprocess.check_output(["./na_3_par_3lvl2_F","-ts_rk_type","5bs","-ts_rtol","1e-8","-ts_atol","1e-8","-n_ens","-1",
                                              "-pulse_type","ARP","-file",fnarr[i],
                                              "-bitstr",init_state[i],
                                              "-b_term",str(b),
                                              "-delta",str(params[0]),
                                              "-pulse_length",str(params[1]),
                                              "-dd_fac",str(ddfac)])
        except:
            pass

        #Read in the QuaC DM
        dm = Qobj(np.loadtxt(fnarr[i]).view(complex),dims=[[2,2,2],[2,2,2]])
        #Remove file
        os.remove(fnarr[i])

        ccz_arp = Qobj([[1,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0],[0,0,0,-1,0,0,0,0],[0,0,0,0,-1,0,0,0],[0,0,0,0,0,-1,0,0],[0,0,0,0,0,0,-1,0],[0,0,0,0,0,0,0,-1]],dims=[[2,2,2],[2,2,2]])
        state_arr=[]
        st=init_state[i]
        for k in range(len(st)):
            if (st[k]=='0'):
                state_arr.append(basis(2,0))
            elif (st[k]=='1'):
                state_arr.append(basis(2,1))
            else:
                print("Invalid input state")
        state=tensor(state_arr[0],state_arr[1],state_arr[2])

        #state=tensor(snot(),snot(),snot())*state
        #Apply phase gates with parameters that we are optimizing
        state = tensor(phasegate(params[2]),qeye(2),qeye(2))*state
        state = tensor(qeye(2),phasegate(params[3]),qeye(2))*state
        state = tensor(qeye(2),qeye(2),phasegate(params[4]))*state
        #Now apply cz_arp
        state = ccz_arp*state

        #Get fidelity wrt quac dm
        fid = fidelity(dm,state)
        fidarr.append(fid)
        earr.append(1-fid)
        print(init_state[i]+' '+str(fid))

    meanf=np.mean(fidarr)
    return 1-meanf
def print_callback(xs):
    print(xs)
print("Optimizing ARP for b = ",str(b))
print("Optimizing Delta, T, and phases")
f.write(str(ddfac)+' ')
#f.write(str(b)+' ')
f.write("Delta_T_phases for b="+str(b)+' ')
default_sp_params = [-0.5,0.2]
default_sp_params = [23,0.54,0,0,0]
res = minimize(fun_sp,default_sp_params,method="nelder-mead",callback=print_callback)

#get the optimal phases
fun_sp(res.x,True)

print("Final Fidelity: ",str(1-res.fun))
f.write(str(1-res.fun)+' ')
print("Final Params: ",str(res.x))
f.write(str(res.x[0])+' ')
f.write(str(res.x[1])+' ')
f.write(str(res.x[2])+' ')
f.write(str(res.x[3])+' ')
f.write(str(res.x[4])+' ')
f.write('\n')
#Final Fidelity:  0.9997463238664505
f.close()

