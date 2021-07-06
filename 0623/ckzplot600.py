from tabulate import tabulate
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
##content=[["qinrui",26,'physics'],["xiaoyu",28,'math']]
##table=tabulate(content,headers=['name','age','major'])
##print(table)
               
hders=["δB","Fidelity","Leakage","φ_1","φ_2","φ_3","φ_4","Δ"]
content=[[1,0.9987,0.9995,0.1896,0.1895,0.1896,0.1896,23.86],
         [2,0.9987,0.9995,0.1894,0.1894,0.1894,0.1894,23.92],
         [5,0.9988,0.9995,0.1894,0.1893,0.1891,0.1894,23.92],
         [10,0.9988,0.9995,0.1888,0.1885,0.1888,0.1887,23.955],
         [20,0.9988,0.9995,0.1879,0.1879,0.1879,0.1879,23.93],
         [30,0.9988,0.9995,0.1871,0.1871,0.1870,0.1871,23.92],
         [40,0.9987,0.9995,0.1864,0.1864,0.1865,0.1864,24.10],
         [50,0.9987,0.9995,0.1858,0.1859,0.1858,0.1860,24.07],
         [60,0.9987,0.9995,0.1853,0.1851,0.1852,0.1854,24.14],
         [100,0.9986,0.9995,0.1843,0.1842,0.1843,0.1844,23.94],
         [150,0.9982,0.9995,0.1857,0.1857,0.1857,0.1854,24.07],
         [200,0.9976,0.9994,0.1900,0.1900,0.1897,0.1898,23.79],
         [300,0.9949,0.9993,0.1896,0.1895,0.1896,0.1896,23.69],]
table=tabulate(content,hders,tablefmt="grid")
print(table)

barr=[val[0] for val in content]
farr=[val[1] for val in content]
print(barr)
print(farr)


plt.plot(barr,farr,'o-',label="Fidelity")
plt.xlabel("δΒ/MHz")
plt.ylabel("Fidelity")
plt.title("F vs δΒ, for B_0=600×2π MHz")
plt.show()
