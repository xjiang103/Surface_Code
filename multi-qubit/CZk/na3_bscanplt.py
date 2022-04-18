import matplotlib.pyplot as plt
import matplotlib
import numpy as np

ax=[10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200]
ay3=[0.9988120300681126,0.9988074155075957,0.9986960976289034,0.9986306845362243,
     0.9985892618225437,0.9985614140061035,0.9985418313335532,0.9985275379809531,
     0.9985168417436887,0.9985085919314034,0.9985021603958261,0.9984969964372215,
     0.9984928289603482,0.9984894350816463,0.9984866087542968,0.9984842283253416,
     0.998482222112535,0.9984805285402724,0.9984790849845574,0.9984778371560283]
plt.plot(ax,ay3)
plt.title("3 atom CZ2,Omega=17MHz")
plt.xlabel("B")
plt.ylabel("Fidelity")
plt.show()
