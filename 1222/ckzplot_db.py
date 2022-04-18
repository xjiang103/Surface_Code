from tabulate import tabulate
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 12})
plt.rcParams["figure.figsize"] = (3.375,2.375)
##content=[["qinrui",26,'physics'],["xiaoyu",28,'math']]
##table=tabulate(content,headers=['name','age','major'])
##print(table)


dbarr2=[197.5,195,190,180,170,150]
dbarr3=[198,196,192,184,176,160]
dbarr4=[198,196,192,184,180,160]
farr3=[0.9992340215542735,0.9992339042183915,0.999227474022634,0.9991943297618133,
       0.9991347620320704,0.9989362532143625]
farr2=[0.9984,0.9984,0.9984,0.9983,0.9982,0.9978]
farr4=[0.9995607107055559,0.9995580052405,0.999549236136917,0.9995398231406177,
       0.9995186761448172,0.99928468405584069]

plt.plot(dbarr2,farr2,'o-',label="2")
plt.plot(dbarr3,farr3,'o-',label="3")
plt.plot(dbarr4,farr4,'o-',label="4")
plt.legend(loc='lower right',  ncol=3, columnspacing=1.5)
plt.ylim([0.9970,1])
plt.xlabel("Min(Β) (MHz)")
plt.ylabel("Fidelity")
#plt.title("F vs min(Β)")

filestr='dbplot.png'
plt.savefig(filestr, bbox_inches='tight',dpi=100)
