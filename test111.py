from qutip import *
from qutip.qip.operations import cnot,snot,phasegate
import numpy as np
from scipy.optimize import minimize
from tabulate import tabulate
import matplotlib.pyplot as plt
import matplotlib

##
##cccz_arp = Qobj(np.diag([1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])
##                   ,dims=[[2,2,2,2],[2,2,2,2]])
####czzz_arp = Qobj(np.diag([1,-1,-1,1,-1,1,1,-1,-1,-1,-1,1,-1,1,1,-1])
####                   ,dims=[[2,2,2,2],[2,2,2,2]])
##    
##    #Get two hadamard state
##state = tensor(snot(),snot(),snot(),snot())*tensor(basis(2,1),basis(2,1),basis(2,1),basis(2,1))
##state = cccz_arp*state
##state = tensor(snot(),qeye(2),qeye(2),qeye(2))*state
##print(state)

state=tensor(basis(2,0),basis(2,0),basis(2,0))
cz=Qobj(np.diag([1,1,1,-1]),dims=[[2,2],[2,2]])
ccz=Qobj(np.diag([1,1,1,1,1,1,1,-1]),dims=[[2,2,2],[2,2,2]])
cccz=Qobj(np.diag([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1]),dims=[[2,2,2,2],[2,2,2,2]])

state=tensor(snot(),snot(),snot())*state
state=tensor(cz,qeye(2))*state
state=tensor(qeye(2),snot(),qeye(2))*state
state=ccz*state
state=tensor(qeye(2),qeye(2),snot())*state
print(state)
