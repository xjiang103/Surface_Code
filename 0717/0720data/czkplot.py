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
            [10, 0.9955443500400134,0.9928250952361852,25.36,0.399],
##            [30,0.9986960976289034,0.9988,26.42,0.463],
##            [50,0.9985892618225437, 0.9995,23.87,0.510],
##            [70,0.9985418313335532,0.9995,23.87,0.510],
##            [90,0.9985168417436887,0.9995,23.87,0.510],
##            [110,0.9985021603958261,0.9995,23.87,0.510],
##            [130,0.9984928289603482,0.9995,23.87,0.510],
##            [150,0.9984866087542968,0.9995,23.87,0.510],
##            [170,0.998482222112535,0.9995,23.87,0.510],
##            [190,0.9984790849845574,0.9995,23.87,0.510],
            [200, 0.9955443500400134,0.9995,23.87,0.510]]
hders3=["k=3,B","Sqrt(F)","Leakage","Δ","T"]
content3=[
            [10, 0.9025624740052018,0.9928250952361852,25.36,0.399],
            [30,0.9164714018825427,0.9988,26.42,0.463],
            [50, 0.9174091939480599,0.9995,23.87,0.510],
            [70, 0.9527130906594035,0.9995,23.87,0.510],
            [90,0.9663602348616578,0.9995,23.87,0.510],
            [110, 0.9716460874596533,0.9995,23.87,0.510],
            [130,0.9781341218002639,0.9995,23.87,0.510],
            [150,0.9825473262017774,0.9995,23.87,0.510],
            [170,0.9856010756756535,0.9995,23.87,0.510],
            [200,0.9886456662928099,0.9995,23.87,0.510]]
hders4=["k=4,B","Sqrt(F)","Leakage","Δ","T"]
content4=[
            [10,0.9340755681006626,0.9928250952361852,25.36,0.399],
            [30, 0.9348953606905765,0.9988,26.42,0.463],
            [50,0.9487388168463682,0.9995,23.87,0.510],
            [70,0.967959198764345,0.9995,23.87,0.510],
            [90,0.9746150708054581,0.9995,23.87,0.510],
            [110, 0.9807308132570051,0.9995,23.87,0.510],
            [130,0.9848111605107612,0.9995,23.87,0.510],
            [150,0.9875396813233368,0.9995,23.87,0.510],
            [170,0.989428163017054,0.9995,23.87,0.510],
            [200, 0.9901622970532064,0.9995,23.87,0.510]]

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
plt.ylim([0.90,1])
plt.xlabel("Β/2π (MHz)")
plt.ylabel("Fidelity")

filestr='fidplot_czk_fix.png'
plt.savefig(filestr, bbox_inches='tight',dpi=100)
plt.show()
