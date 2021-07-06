from qutip import *
from qutip.qip.operations import cnot,snot
import numpy as np

zero=basis(2,0)
one=basis(2,1)

state=one

one_one=one*one.dag()

print(one_one)

H=-1j*one_one

times=np.linspace(0,100,10000)

result=sesolve(H,state,times)

res=result.states[-1]

print(res)


