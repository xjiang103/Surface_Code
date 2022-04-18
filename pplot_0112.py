import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy.optimize import curve_fit

def func(x,a,b):
    return a*x+b

matplotlib.rcParams.update({'font.size': 16})
f=open("phase_0112.txt","r")
xa=[]
y1a=[]
y2a=[]
y3a=[]
y2appi=[]
y2ampi=[]
flag1=0
flag2=0
for i in range(100):
    r=f.readline()
    x,y1,y2,y3=r.split()
    print(x)
    print(y1)
    xa.append(float(x))
    y1a.append(float(y1))
    y2a.append(float(y2))
    y3a.append(float(y3))
    y2appi.append(float(y2)+np.pi)
    y2ampi.append(float(y2)-np.pi)
    if((float(y1)>float(y2)+np.pi)&(flag1==0)):
        print("x1=")
        print(float(x))
        flag1=1
    if((float(y3)<float(y2)-np.pi)&(flag2==0)):
        print("x2=")
        print(float(x))
        flag2=1
zx=xa[0:49]
z1=y1a[0:49]
z2=y2a[0:49]
z3=y3a[0:49]

popt, pcov = curve_fit(func,zx,z1)
print("parameter for 1111")
print(popt)
popt, pcov = curve_fit(func,zx,z2)
print("parameter for 1110")
print(popt)
popt, pcov = curve_fit(func,zx,z3)
print("parameter for 1100")
print(popt)

#1111 and 1110 : x=0.220747
#1110 and 1100 : x=0.223737

#try 0.1 data-data coupling

#try npi difference

plt.plot(xa,y1a,'-',label="1111")
plt.plot(xa,y2a,'-',label="1110")
plt.plot(xa,y3a,'-',label="1100")
plt.plot(xa,y2appi,'--',label="1110+pi")
plt.plot(xa,y2ampi,'--',label="1110-pi")
plt.xlabel("δτ")
plt.ylabel("phase")

plt.legend()
plt.savefig("phaseplot1_0117.png")
plt.show()

