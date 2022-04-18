from tabulate import tabulate
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 20})
##content=[["qinrui",26,'physics'],["xiaoyu",28,'math']]
##table=tabulate(content,headers=['name','age','major'])
##print(table)


dbarr1=[0,2,5,10,20,40,60,90]
farr1=[0.9984,0.9984,0.9984,0.9984,0.9983,0.9980,0.9975,0.9942]

plt.plot(dbarr1,farr1,'o-')
#plt.legend(loc='right')
plt.xlabel("dΒ/MHz")
plt.ylabel("Fidelity")
plt.title("F vs dΒ, 3atom, for B0=200,B_list=[B0-dB,B0,B0+dB]")
plt.show()
