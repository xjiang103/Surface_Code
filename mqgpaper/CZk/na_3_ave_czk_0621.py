from scipy.optimize import minimize
import subprocess
import numpy as np
import uuid
import os
from qutip import *
from qutip.qip.operations import snot,phasegate
import argparse
from joblib import Parallel, delayed

import random


parser = argparse.ArgumentParser(description='Optimize Pulses')
parser.add_argument('-b','--b_couple', help='Rydberg Coupling Strength', required=True, type=float)
parser.add_argument('-pt','--pulse_type', help='pulse type', required=True, type=int)
args = vars(parser.parse_args())
b = args["b_couple"]
ptnum = args["pulse_type"]
#ist=args["initial_state"]
#print("initial state is "+str(ist))
ddfac=0
typestr="cccz"
default_sp_params = [0.5045,0.222222]
ptype="SP"
para1str="deltat"

if (ptnum==1):
    ptype="ARP"
    randa=random.uniform(18,28)
    randb=random.uniform(0.2,0.8)
    default_sp_params = [randa,randb]
    para1str="-pulse_length"
elif (ptnum==2):
    ptype="SP"
    para1str="-deltat"
    randa=random.uniform(0.2,0.8)
    randb=random.uniform(0.1,0.4)
    default_sp_params = [randa,randb]
else:
    print("pulse type error, must be 1 or 2")

unique_file = str(uuid.uuid4())[0:8]
file_name = "dm_"+unique_file+".dat" #Allow us to run in parallel
params = [23]
#b=100
f=open("arp_3_F_ave.txt","a")
f.write('\n')
f.write(ptype+' ')
f.write(str(default_sp_params[0])+' ')
f.write(str(default_sp_params[1])+' ')
#f.write("initial state="+str(ist)+'\n')
f.write(typestr+' ')

#--------------------------------------------------
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

#add the additional superposition state
for j in range(8):
    ave_state+=(1/np.sqrt(8))*init_arr[j]
init_arr.append(ave_state)
str_arr.append("xxx")

#--------------------------------------------------------

#init_state=init_arr[ist]

#print("init_state is the "+str_arr[ist]+"th basis state. ")
#print("init state:")
#print(init_state)
#----------------------------------------------

def run_job(i,params):
        unique_file = str(uuid.uuid4())[0:8]
        file_name = "dm_"+unique_file+".dat" #Allow us to run in parallel
        #Run QuaC
        try:
            output = subprocess.check_output(["./na_3_F_par","-ts_rk_type","5bs","-ts_rtol","1e-8","-ts_atol","1e-8","-n_ens","-1",
                                              "-pulse_type",ptype,"-file",file_name,
                                              "-bitstr",str_arr[i],
                                              "-b_term",str(b),
                                              "-delta",str(params[0]),
                                              para1str,str(params[1]),
                                              "-dd_fac",str(ddfac)])
        except:
            output = subprocess.check_output(["./na_3_F_par","-ts_rk_type","5bs","-ts_rtol","1e-8","-ts_atol","1e-8","-n_ens","-1",
                                              "-pulse_type",ptype,"-file",file_name,
                                              "-bitstr",str_arr[i],
                                              "-b_term",str(b),
                                              "-delta",str(params[0]),
                                              para1str,str(params[1]),
                                              "-dd_fac",str(ddfac)])

            print("fail")
            pass

        #Read in the QuaC DM
        #dm_arr.append(Qobj(np.loadtxt(file_name).view(complex),dims=[[2,2,2],[2,2,2]]))

        dm = Qobj(np.loadtxt(file_name).view(complex),dims=[[2,2,2],[2,2,2]])
        #Remove file
        os.remove(file_name)
        return dm

def fun_sp(params,final_run=None):

    results = Parallel(n_jobs=9,backend="loky")(delayed(run_job)(i,params) for i in range(9))
    #QUTIP to get perfect circuit
    res = minimize(qutip_phase,[0,0,0],method="COBYLA",args=(results))

    fid = 1-res.fun
    print(fid)
    if(final_run):
        print("Phase: ",res.x)
        f.write(str(res.x[0])+' ')
        f.write(str(res.x[1])+' ')
        f.write(str(res.x[2])+' ')
    return 1-fid
def print_callback(xs):
    print(xs)

def qutip_phase(params,dms):
    #define cz_arp and czz arp
    ccz_arp = Qobj([[1,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0],[0,0,0,-1,0,0,0,0],[0,0,0,0,-1,0,0,0],[0,0,0,0,0,-1,0,0],[0,0,0,0,0,0,-1,0],[0,0,0,0,0,0,0,-1]],dims=[[2,2,2],[2,2,2]])
    czz_arp = Qobj([[1,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1]],dims=[[2,2,2],[2,2,2]])
    fid=0
    fid_m=1
    #leakage=0
    fid_tmp=0
    f.write('\n')
    for i in range(9):
        state=init_arr[i]
        #Apply phase gates with parameters that we are optimizing
        state = tensor(phasegate(params[0]),qeye(2),qeye(2))*state
        state = tensor(qeye(2),phasegate(params[1]),qeye(2))*state
        state = tensor(qeye(2),qeye(2),phasegate(params[2]))*state

        #Now apply cz_arp
        state = czz_arp*state

        #Get fidelity wrt quac dm
        fid_tmp = (fidelity(dms[i],state))
        #leak_tmp = (np.trace(dms[i])).real
        #print(str(i)+' '+str(fid_tmp))
        #print(str(i)+' '+str(leak_tmp))
        fid=fid+fid_tmp/9.0
        fid_m=fid_m*fid_tmp
        #leakage=leakage+leak_tmp/17.0
        #f.write(str(fid_tmp)+'\n')
    fid_m=fid_m/fid_tmp
    lambda1=1-(1-fid_m)/(1-fid_tmp*fid_m)
    fg=1/9+8/9*fid_m*fid_tmp
    f_final=lambda1*fg+fid*(1-lambda1)
    #print("lambda is "+str(lambda1)+", F="+str(f_final)+"\n")
    #f.write('\n')
    return 1-f_final

def fun_arp(delta):
    #NOT COMPLETED!
    return 1-fid


print("Optimizing ARP for b = ",str(b))
print("Optimizing Delta, T, and phases")
f.write(str(ddfac)+' ')
#f.write(str(b)+' ')
f.write("Delta_T_phases for b="+str(b)+' ')

res = minimize(fun_sp,default_sp_params,method="nelder-mead",callback=print_callback)

#get the optimal phases
fun_sp(res.x,True)

print("Final Fidelity: ",str(1-res.fun))
f.write(str(1-res.fun)+' ')
print("Final Params: ",str(res.x))
f.write(str(res.x[0])+' ')
f.write(str(res.x[1])+' ')
f.write('\n')
#Final Fidelity:  0.9997463238664505
f.close()

