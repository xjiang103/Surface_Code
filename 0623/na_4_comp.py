from tabulate import tabulate
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
##content=[["qinrui",26,'physics'],["xiaoyu",28,'math']]
##table=tabulate(content,headers=['name','age','major'])
##print(table)
hders=["B","F_QuaC","F_QuTiP","Leak_QuaC","Leak_QuTiP","φ_1","φ_2","φ_3","φ_4","Δ"]
content=[[200,0.9907,0.9891,0.9995,0.9994,0.5726,0.5726,0.5728,0.5725,23.86],
         [300,0.9957,0.9945,0.9995,0.9994,0.3805,0.3805,0.3803,0.3805,23.94],
         [400,0.9975,0.9964,0.9995,0.9994,0.2849,0.2848,0.2849,0.2850,24.02],
         [500,0.9983,0.9973,0.9995,0.9994,0.2277,0.2276,0.2277,0.2276,23.86],
         [600,0.9987,0.9978,0.9995,0.9994,0.1896,0.1895,0.1896,0.1896,23.86]]

table=tabulate(content,hders,tablefmt="grid")
print(table)

barr=[val[0] for val in content]
fquacarr=[val[1] for val in content]
fqutiparr=[val[2] for val in content]


plt.plot(barr,fquacarr,'o-',label="QuaC")
plt.plot(barr,fqutiparr,'o-',label="QuTiP")
plt.legend(loc='right')
plt.xlabel("B/MHz")
plt.ylabel("Fidelity")
plt.title("F vs Β, with 4 qubits")
plt.show()
