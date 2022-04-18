import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
f=open("phase_0105.txt","r")
xa=[]
y1a=[]
y2a=[]
y3a=[]
for i in range(18):
    r=f.readline()
    x,y1,y2=r.split()
    print(x)
    print(y1)
    xa.append(float(x))
    y1a.append(float(y1))
    y2a.append(float(y2))
    y3a.append(float(y1)-float(y2))


plt.plot(xa,y1a,'o-',label="1111")
plt.plot(xa,y2a,'o-',label="1110")
plt.plot(xa,y3a,'o-',label="1111-1110")

plt.legend()

plt.show()

