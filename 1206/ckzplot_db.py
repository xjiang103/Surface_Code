from tabulate import tabulate
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 20})
##content=[["qinrui",26,'physics'],["xiaoyu",28,'math']]
##table=tabulate(content,headers=['name','age','major'])
##print(table)


dbarr1=[0,2,5,10,20,40,60,90]
dbarr4=[0,2,5,10,20,40]
farr3=[0.9992,0.9992,0.9992,0.9991,0.9989,0.9978,0.9934,0.9684]
farr2=[0.9984,0.9984,0.9984,0.9984,0.9983,0.9980,0.9975,0.9942]
farr4=[0.9991,0.9991,0.9990,0.9989,0.9980,0.9305]

plt.plot(dbarr1,farr2,'o-',label="3 atom")
plt.plot(dbarr1,farr3,'o-',label="4 atom")
plt.plot(dbarr4,farr4,'o-',label="5 atom")
plt.legend(loc='right')
plt.xlabel("dΒ/MHz")
plt.ylabel("Fidelity")
plt.title("F vs dΒ")
plt.show()
