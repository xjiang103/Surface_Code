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
            [10,0.9988120300681126,0.9928250952361852,25.36,0.399],
            [30,0.9986960976289034,0.9988,26.42,0.463],
            [50,0.9985892618225437, 0.9995,23.87,0.510],
            [70,0.9985418313335532,0.9995,23.87,0.510],
            [90,0.9985168417436887,0.9995,23.87,0.510],
            [110,0.9985021603958261,0.9995,23.87,0.510],
            [130,0.9984928289603482,0.9995,23.87,0.510],
            [150,0.9984866087542968,0.9995,23.87,0.510],
            [170,0.998482222112535,0.9995,23.87,0.510],
            [190,0.9984790849845574,0.9995,23.87,0.510],
            [200,0.9984778371560283,0.9995,23.87,0.510]]
hders3=["k=3,B","Sqrt(F)","Leakage","Δ","T"]
content3=[
            [10,0.9928250948202235,0.9928250952361852,25.36,0.399],
            [30,0.9932499185195285,0.9988,26.42,0.463],
            [50,0.9934244641146552,0.9995,23.87,0.510],
            [70,0.9935145222556008,0.9995,23.87,0.510],
            [90,0.9935696950636718,0.9995,23.87,0.510],
            [110,0.9936073517877018,0.9995,23.87,0.510],
            [130,0.9936508079531619,0.9995,23.87,0.510],
            [150,0.9936767626764029,0.9995,23.87,0.510],
            [170,0.993696643047185,0.9995,23.87,0.510],
            [190,0.9937123329195241,0.9995,23.87,0.510],
            [200,0.9937189943843413,0.9995,23.87,0.510]]
hders4=["k=4,B","Sqrt(F)","Leakage","Δ","T"]
content4=[
            [10,0.9967447753389103,0.9928250952361852,25.36,0.399],
            [30, 0.9950215776565405,0.9988,26.42,0.463],
            [50,0.996620785404963,0.9995,23.87,0.510],
            [70,0.9974636221088323,0.9995,23.87,0.510],
            [90,0.9969420567053532,0.9995,23.87,0.510],
            [110,0.9972325300890108,0.9995,23.87,0.510],
            [130,0.9974070149530665,0.9995,23.87,0.510],
            [150,0.9970790476733762,0.9995,23.87,0.510],
            [170,0.9973872120875011,0.9995,23.87,0.510],
            [190,0.9962999165760821,0.9995,23.87,0.510]]

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
plt.ylim([0.99,1])
plt.xlabel("Β/2π (MHz)")
plt.ylabel("Fidelity")

filestr='fidplot_czk.png'
plt.savefig(filestr, bbox_inches='tight',dpi=100)
plt.show()
