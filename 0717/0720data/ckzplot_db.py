from tabulate import tabulate
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 12})
plt.rcParams["figure.figsize"] = (3.375,2.375)
##content=[["qinrui",26,'physics'],["xiaoyu",28,'math']]
##table=tabulate(content,headers=['name','age','major'])
##print(table)


dbarr2=[200,198,195,190,185,175,150]
dbarr3=[200,198,195,190,185,175,150]
dbarr4=[200,196,190,180,170,160,150]
farr2=[0.991432742737744,
       0.9914275488285758,
       0.9914083482025744,
       0.9913366899276851,
       0.9912165653781774,
       0.9908181027809533,
       0.9886416208168587]

farr3=[0.9929650192943118,
       0.9929792098877821,
       0.9929768218304468,
       0.9928946697790965,
       0.9927184951668337,
       0.9921041852285526,
       0.9889692645554706]

farr4=[0.9954379345528093,
       0.9953822324621768,
       0.9952374089352791,
       0.9946269824977482,
       0.993583715037088,
       0.9920748104833814,
       0.9900388967179834]

plt.plot(dbarr2,farr2,'o-',label="2")
plt.plot(dbarr3,farr3,'o-',label="3")
plt.plot(dbarr4,farr4,'o-',label="4")
plt.legend(loc='lower right',  ncol=3, columnspacing=1.5)
plt.ylim([0.9780,0.996])
plt.xlabel("Min(Β/2π) (MHz)")
plt.ylabel("Fidelity")
#plt.title("F vs min(Β)")

filestr='dbplot_ckz_fix.png'
plt.savefig(filestr, bbox_inches='tight',dpi=100)
