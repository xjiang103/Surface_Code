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
parser.add_argument('-deltab','--deltab_couple', help='Variation in Rydberg Coupling Strength', required=True, type=float)
args = vars(parser.parse_args())
b = args["b_couple"]
deltab=args["deltab_couple"]

ddfac=1
typestr="cccz"


unique_file = str(uuid.uuid4())[0:8]
file_name = "dm_"+unique_file+".dat" #Allow us to run in parallel
params = [23]
#b=100
f=open("arp_4_F_ave_db.txt","a")
f.write('\n')
#f.write("initial state="+str(ist)+'\n')
f.write(typestr+' ')

params_pulse=[23.95157083,0.53117221]
params_phase=[0.56335345,0.5633543,0.56319147,0.56335745]
#--------------------------------------------------
def ia(init_state):
    state_arr=[]
    for k in range(len(init_state)):
        if (init_state[k]=='0'):
            state_arr.append(basis(2,0))
        elif (init_state[k]=='1'):
            state_arr.append(basis(2,1))
    state=tensor(state_arr[0],state_arr[1],state_arr[2],state_arr[3])
    return state
init_arr=[]
str_arr=[]

#init_arr is an arry of initial states
#str_arr is an arry of strings, to be used in calling the QuaC program

#enumerate the basis states
for j1 in range(2):
    for j2 in range(2):
        for j3 in range(2):
            for j4 in range(2):
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
                if (j4==0):
                    initstr=initstr+'0'
                elif (j4==1):
                    initstr=initstr+'1'
                str_arr.append(initstr)
                init_arr.append(ia(initstr))
ave_state=0

#add the additional superposition state
for j in range(16):
    ave_state+=(1/np.sqrt(16))*init_arr[j]
init_arr.append(ave_state)
str_arr.append("xxxx")
#print("len="+str(init_arr))

#--------------------------------------------------------

#init_state=init_arr[ist]

#print("init_state is the "+str_arr[ist]+"th basis state. ")
#print("init state:")
#print(init_state)
#----------------------------------------------
params=params_pulse
def run_job(i,params):
        unique_file = str(uuid.uuid4())[0:8]
        file_name = "dm_"+unique_file+".dat" #Allow us to run in parallel
        #Run QuaC
        try:
            output = subprocess.check_output(["./na_4_F_par_db","-ts_rk_type","5bs","-ts_rtol","1e-8","-ts_atol","1e-8","-n_ens","-1",
                                              "-pulse_type","ARP","-file",file_name,
                                              "-bitstr",str_arr[i],
                                              "-b_term",str(b),
                                              "-delta_b_term",str(deltab),
                                              "-delta",str(params[0]),
                                              "-pulse_length",str(params[1]),
                                              "-dd_fac",str(ddfac)])
        except:
            output = subprocess.check_output(["./na_4_F_par_db","-ts_rk_type","5bs","-ts_rtol","1e-8","-ts_atol","1e-8","-n_ens","-1",
                                              "-pulse_type","ARP","-file",file_name,
                                              "-bitstr",str_arr[i],
                                              "-b_term",str(b),
                                              "-delta_b_term",str(deltab),
                                              "-delta",str(params[0]),
                                              "-pulse_length",str(params[1]),
                                              "-dd_fac",str(ddfac)])

            print("fail")
            pass

        #Read in the QuaC DM
        #dm_arr.append(Qobj(np.loadtxt(file_name).view(complex),dims=[[2,2,2],[2,2,2]]))

        dm = Qobj(np.loadtxt(file_name).view(complex),dims=[[2,2,2,2],[2,2,2,2]])
        #Remove file
        os.remove(file_name)
        return dm

results = Parallel(n_jobs=17,backend="loky")(delayed(run_job)(i,params) for i in range(17))
dms=results
#QUTIP to get perfect circuit

cccz_arp = Qobj(np.diag([1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])
               ,dims=[[2,2,2,2],[2,2,2,2]])
czzz_arp = Qobj(np.diag([1,-1,-1,1,-1,1,1,-1,-1,-1,-1,1,-1,1,1,-1])
               ,dims=[[2,2,2,2],[2,2,2,2]])
fid=0

params=params_phase
for i in range(17):
    state=init_arr[i]
    #Apply phase gates with parameters that we are optimizing
    state = tensor(phasegate(params[0]),qeye(2),qeye(2),qeye(2))*state
    state = tensor(qeye(2),phasegate(params[1]),qeye(2),qeye(2))*state
    state = tensor(qeye(2),qeye(2),phasegate(params[2]),qeye(2))*state
    state = tensor(qeye(2),qeye(2),qeye(2),phasegate(params[3]))*state

    #Now apply cz_arp
    state = cccz_arp*state

    #Get fidelity wrt quac dm
    fid_tmp = fidelity(dms[i],state)
    fid=fid+fid_tmp/17.0

print("Optimizing ARP for b = ",str(b))
print("Delta_b = ",str(deltab))
f.write(str(b)+' '+str(deltab)+' ')
#f.write(str(b)+' ')
f.write("Delta_T_phases for b="+str(b)+' ')
f.write(str(params_pulse[0])+' '+str(params_pulse[1])+' ')
f.write(str(params_phase[0])+' '+str(params_phase[1])+' '+str(params_phase[2])+' '+str(params_phase[3])+' ')

print("Final Fidelity: "+str(fid))
f.write(str(fid))
f.write('\n')
#Final Fidelity:  0.9997463238664505
f.close()

