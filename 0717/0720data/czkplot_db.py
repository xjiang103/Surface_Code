from tabulate import tabulate
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib
matplotlib.rcParams.update({'font.size': 12})
plt.rcParams["figure.figsize"] = (3.375,2.375)
##content=[["qinrui",26,'physics'],["xiaoyu",28,'math']]
##table=tabulate(content,headers=['name','age','major'])
##print(table)


dbarr=[0e-5,1e-5,2e-5,4e-5,6e-5,8e-5,10e-5]
farr2=[0.9955443515741711,
       0.9954863591312847,
       0.995298495828996,
       0.9945323834798147,
       0.9932476759704801,
       0.9914464530840894,
       0.9891313910182727]
farr3=[0.9886456849279844,
       0.9885647735453649,
       0.9884803240503013,
       0.9882900223595863,
       0.9880892329962019,
       0.9878615843815082,
       0.9876122113577601]
farr4=[0.990162310397958,
       0.990044273580095,
       0.9899343107182296,
       0.9896280685010069,
       0.9892446268243918,
       0.9887608428880799,
       0.988246348844815]
plt.ticklabel_format(style='sci', axis='x')
plt.plot(dbarr,farr2,'o-',label="2")
plt.plot(dbarr,farr3,'o-',label="3")
plt.plot(dbarr,farr4,'o-',label="4")
plt.legend(loc='lower right',  ncol=3, columnspacing=1.5)
plt.ylim([0.9830,1])
plt.xlabel("T-T/C-T coupling")
plt.ylabel("Fidelity")
#plt.title("F vs min(Β)")

class ScalarFormatterClass(ScalarFormatter):
   def _set_format(self):
      self.format = "%1.2f"
ax = plt.gca()
xScalarFormatter = ScalarFormatterClass(useMathText=True)
xScalarFormatter.set_powerlimits((0,0))
ax.xaxis.set_major_formatter(xScalarFormatter)

filestr='czk_dbplot_fix.png'
plt.savefig(filestr, bbox_inches='tight',dpi=100)
plt.show()
