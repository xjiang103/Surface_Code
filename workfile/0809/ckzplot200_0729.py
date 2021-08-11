from tabulate import tabulate
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
##content=[["qinrui",26,'physics'],["xiaoyu",28,'math']]
##table=tabulate(content,headers=['name','age','major'])
##print(table)
n0=4

hders2=["k=2,B","Sqrt(F)","Leakage","Δ","T"]
content2=[[100,0.9722,0.9891,28.46,0.346],
            [200,0.9915,0.9973,25.36,0.399],
            [300,0.9956,0.9979,25.02,0.408],
            [400,0.9970,0.9988,26.42,0.463],
            [500,0.9979,0.9992,25.49,0.486],
            [600,0.9984,0.9995,23.87,0.510]]
hders3=["k=3,B","Sqrt(F)","Leakage","Δ","T"]
content3=[[100,0.9892,0.9910,29.54,0.345],
            [200,0.9935,0.9973,27.12,0.431],
            [300,0.9966,0.9985,27.35,0.440],
            [400,0.9979,0.9989,27.26,0.460],
            [500,0.9984,0.9993,24.79,0.512],
            [600,0.9988,0.9995,24.03,0.525]]
hders4=["k=4,B","Sqrt(F)","Leakage","Δ","T"]
content4=[[100,0.9877,0.994,30.30,0.356],
            [200,0.9959,0.9984,28.66,0.427]]

barr=[1,2,3,4]
farr=[0.9955,0.9956,0.9966,0.9978]
plt.plot(barr,farr,'o-',label="B=200")
plt.legend(loc='right')
plt.xlabel("k")
plt.ylabel("Fidelity")
plt.title("F vs k, for B=200")
plt.show()

if(n0==2):
    content=content2
    hders=hders2
if(n0==3):
    content=content3
    hders=hders3
if(n0==4):
    content=content4
    hders=hders4
table=tabulate(content,hders,tablefmt="grid")
print(table)

barr=[val[0] for val in content]
farr=[val[1] for val in content]
leakarr=[val[2] for val in content]
print(barr)
print(farr)


##plt.plot(barr,farr,'o-',label="Fidelity")
##plt.plot(barr,leakarr,'o-',label="Leakage")
##plt.legend(loc='right')
##plt.xlabel("Β/MHz")
##plt.ylabel("Fidelity")
##plt.title("F vs Β, for k="+str(n0))
##plt.show()
