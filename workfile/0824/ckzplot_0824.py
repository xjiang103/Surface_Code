from tabulate import tabulate
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
##content=[["qinrui",26,'physics'],["xiaoyu",28,'math']]
##table=tabulate(content,headers=['name','age','major'])
##print(table)


dbarr1=[0.25,1.25,2.5,5,10,15,20]
farr1=[0.9959,0.9959,0.9958,0.9950,0.9895,0.9812,0.96]

plt.plot(dbarr1,farr1,'o-',label="Use the same parameter")
plt.legend(loc='right')
plt.xlabel("dΒ/MHz")
plt.ylabel("Fidelity")
plt.title("F vs dΒ, 5atom, for B0=200,B_list=[B0-4dB,B0-3dB,...,B0+5dB]")
plt.show()
