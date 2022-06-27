import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
f=open("f3swp1_1113.txt","r")
xa=[]
ya=[]
for i in range(80):
    r=f.readline()
    x,y,sp,sn=r.split()
    print(str(i)+" "+x+" "+y)
    xa.append(float(x))
    ya.append(float(y))


plt.plot(xa,ya,'o-')
plt.xlabel("fcenter/(Ω0/2π))")
plt.ylabel('Error')
plt.title("2-photon, t=1π/Ω")
plt.show()

