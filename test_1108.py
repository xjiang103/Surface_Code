from scipy.optimize import minimize
import subprocess
import numpy as np
import uuid
import os
from qutip import *
from qutip.qip.operations import snot,phasegate
import argparse

one=basis(3,1)
zero=basis(3,0)
two=basis(3,2)

state=(1/np.sqrt(3))*(one+zero+two)
print(state)

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
            print(initstr)
            init_arr.append(ia(initstr))
ave_state=0
for j in range(8):
    ave_state+=(1/np.sqrt(8))*init_arr[j]
init_arr.append(ave_state)
print(init_arr[8])
            
