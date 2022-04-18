from scipy.optimize import minimize
import subprocess
import numpy as np
import uuid
import os
from qutip import *
from qutip.qip.operations import snot,phasegate
import argparse

dtarr=[]
arr1=[]
arr2=[]
arr3=[]
ar1tmp=0
ct1tmp=0
ar2tmp=0
ct2tmp=0
ar3tmp=0
ct3tmp=0


f1=open("phase_0112.txt","w")
for kkk in range(100):
    dt_in=0.15+kkk*0.004
    params=0.5
    
    subprocess.run(["./na_0105",
                                    "-deltat",str(dt_in),
                                    "-bitstr","1111"])
    data = np.loadtxt("na_0105.txt")
    dtarr.append(dt_in)
    
##    ar1tmp=0
##    ct1tmp=0

    datatmp=data[:,9][-1]
    datatmp=datatmp+np.pi
    if((datatmp+ct1tmp*2*np.pi)>ar1tmp):
        ar1tmp=datatmp+ct1tmp*2*np.pi
        arr1.append(ar1tmp)
    else:
        ct1tmp=ct1tmp+1
        ar1tmp=datatmp+2*np.pi*ct1tmp
        arr1.append(ar1tmp)
        

    output = subprocess.run(["./na_0105",
                                    "-deltat",str(dt_in),
                                    "-bitstr","1110"])
    data = np.loadtxt("na_0105.txt")

##    ar2tmp=0
##    ct2tmp=0
    
    datatmp=data[:,10][-1]
    datatmp=datatmp+np.pi
    if(datatmp+ct2tmp*2*np.pi>ar2tmp):
        ar2tmp=datatmp+ct2tmp*2*np.pi
        arr2.append(ar2tmp)
    else:
        ct2tmp=ct2tmp+1
        ar2tmp=datatmp+2*np.pi*ct2tmp
        arr2.append(ar2tmp)

    output = subprocess.run(["./na_0105",
                                    "-deltat",str(dt_in),
                                    "-bitstr","1100"])
    data = np.loadtxt("na_0105.txt")

##    ar3tmp=0
##    ct3tmp=0
    
    datatmp=data[:,11][-1]
    datatmp=datatmp+np.pi
    if(datatmp+ct3tmp*2*np.pi>ar3tmp):
        ar3tmp=datatmp+ct3tmp*2*np.pi
        arr3.append(ar3tmp)
    else:
        ct3tmp=ct3tmp+1
        ar3tmp=datatmp+2*np.pi*ct3tmp
        arr3.append(ar3tmp)


    f1.write(str(dtarr[kkk])+' '+str(arr1[kkk])+' '+str(arr2[kkk])+' '+str(arr3[kkk])+'\n')

f1.close()

    
    
    
