from tabulate import tabulate
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 12})
plt.rcParams["figure.figsize"] = (3.375,2.375)
##content=[["qinrui",26,'physics'],["xiaoyu",28,'math']]
##table=tabulate(content,headers=['name','age','major'])
##print(table)
n0=4

hders2=["k=2,B","Sqrt(F)","Leakage","Δ","T"]
content2=[
            [200,0.9914,0.9973,25.36,0.399],
            [300,0.9950,0.9979,25.02,0.408],
            [400,0.9969,0.9988,26.42,0.463],
            [500,0.9979,0.9992,25.49,0.486],
            [600,0.9983,0.9995,23.87,0.510]]
hders3=["k=3,B","Sqrt(F)","Leakage","Δ","T"]
content3=[
            [200,0.9929650264903672,0.9973,27.12,0.431],
            [300,0.9961982138781537,0.9985,27.35,0.440],
            [400,0.997487461941445,0.9989,27.26,0.460],
            [500,0.998211035756138,0.9993,24.79,0.512],
            [600,0.9986107329315749,0.9995,24.03,0.525]]
hders4=["k=4,B","Sqrt(F)","Leakage","Δ","T"]
content4=[
            [200,0.9954379281868172,0.9973,27.12,0.431],
            [300,0.9974215545538365,0.9985,27.35,0.440],
            [400,0.9981945268300926,0.9989,27.26,0.460],
            [500,0.9984945268300926,0.9993,24.79,0.512],
            [600,0.9987945268300926,0.9995,24.03,0.525]]

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

barr=[val[0] for val in content2]
farr=[val[1] for val in content2]
#leakarr=[val[2] for val in content]
plt.plot(barr,farr,'o-',label="2")

barr=[val[0] for val in content3]
farr=[val[1] for val in content3]
#leakarr=[val[2] for val in content]
plt.plot(barr,farr,'o-',label="3")

barr=[val[0] for val in content4]
farr=[val[1] for val in content4]
#leakarr=[val[2] for val in content]
plt.plot(barr,farr,'o-',label="4")
#plt.plot(barr,leakarr,'o-',label="Leakage")

plt.legend(loc='lower right', prop={'size': 8}, ncol=3)
plt.xlabel("Β/2π (MHz)")
plt.ylabel("Fidelity")

filestr='fidplot_fix.png'
plt.savefig(filestr, bbox_inches='tight',dpi=100)
plt.show()
