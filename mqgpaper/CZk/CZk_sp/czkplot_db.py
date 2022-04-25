from tabulate import tabulate
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 12})
plt.rcParams["figure.figsize"] = (3.375,2.375)
##content=[["qinrui",26,'physics'],["xiaoyu",28,'math']]
##table=tabulate(content,headers=['name','age','major'])
##print(table)


dbarr=[0,0.05,0.1,0.15,0.2,0.25]
farr2=[0.9988230567464242,0.9986701384551281,0.9975164959274561,0.9955716490752029,
       0.9930156664737482,0.9899933136907312]
farr3=[0.9928250952361852 ,0.9925117275367986,0.9920964528061972,
       0.9916012908094677,0.991044624724391 ,0.9898036291091025]
farr4=[0.998,0.9975,0.9965,0.995,0.9925,0.990]

plt.plot(dbarr,farr2,'o-',label="2")
plt.plot(dbarr,farr3,'o-',label="3")
plt.plot(dbarr,farr4,'o-',label="4")
plt.legend(loc='lower right',  ncol=3, columnspacing=1.5)
plt.ylim([0.9850,1])
plt.xlabel("T-T coupling/C-T coupling")
plt.ylabel("Fidelity")
#plt.title("F vs min(Β)")

filestr='czk_dbplot.png'
plt.savefig(filestr, bbox_inches='tight',dpi=100)
