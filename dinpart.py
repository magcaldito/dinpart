import matplotlib.pyplot as plt
import numpy as np
import random

def plotvector(v):
    x=np.array((0,v[0]))
    y=np.array((0,v[1]))
#    print(x,y)
    plt.plot(x,y)
    
def plotvectorpunto(v,punto):
    x=np.array((punto[0],punto[0]+v[0]))
    y=np.array((punto[1],punto[1]+v[1]))
#    print(x,y)
    plt.plot(x,y)

def modulo(v):
    modulo=(v[0]**2+v[1]**2)**0.5
    return modulo

# a=np.array((2,5))
# b=np.array((3,2))

# print (a,b,a+b,a-b)

# fig = plt.figure(figsize=(8,8))
# ax = fig.add_subplot(1, 1, 1)
# ax.spines['left'].set_position(('data',0))
# ax.spines['bottom'].set_position(('data',0))
# ax.spines['right'].set_color('none')
# ax.spines['top'].set_color('none')
# ax.xaxis.set_ticks_position('bottom')
# ax.yaxis.set_ticks_position('left')
# plt.xlim(-10,10)
# plt.ylim(-10,10)

# plotvector(a)
# plotvectorpunto(a,b)
# plotvectorpunto(a-b,b)
# #plotvectorpunto(b,a)

# #plotvector(b)
# #plotvector(a+b)
# #plotvector(a-b)

# plt.show()

#t1=(2,3)
#t2=(1,1)
#l1=[2,3]
#l2=[1,1]
#print(t1,t2,t1+t2)
#print(l1,l2,l1+l2)

G= 6.674E-11            # Gravity constant 
ME=5.97219E24           # Mass of the earth
MS=1.9891E30            # Mass of the sun
MM=7.3476E22            # Mass of the moon
DES= 1.50E11            # Distance to the sun
DEM=3.844E8             # Distance to the moon
VES= 29784.8            # Velocity earth around sun
VME=1023                # Velocity moon around earth

escala=0.0000001


n=8                     # no. particles
nt=500                  # no. time steps
ts=5000                 # time step
t=np.zeros(nt)          # time vector

m=np.zeros(n)           # particle masses vector  


d=np.zeros((n,nt,2))    # position vectors 
v=np.zeros((n,nt,2))    # velocity vectors
a=np.zeros((n,nt,2))    # acceleration vectors
f=np.zeros((n,nt,2))    # force vectors
vaux=np.zeros((n,nt,2))    # force vectors

print('---------------------------------')

#print(modulo([1,1]))

for i in range (nt):
    t[i]=i*ts

for i in range(n):          # Condciones iniciales aleatorias
    m[i]=random.uniform(8000,12000)*ME
    d[i][0]=[random.uniform(0.1,0.2)*DES,random.uniform(0.1,0.2)*DES]
    v[i][0]=[random.uniform(-1,1)*VES,random.uniform(-1,1)*VES]
    print(d[i][0][1])

# m[0]=MS
# m[1]=ME
# m[2]=MM

# d[0,0]=[0,0]
# d[1,0]=[0,DES]
# d[2,0]=[0,DES+DEM]

# v[0,0]=[0,0]
# v[1,0]=[VES,0]
# v[2,0]=[VES+VME,0]

# for k in range(nt):
#    for i in range(n):
#        f[i][k]=[0,0]
#        for j in range(n):
#            if j==i:
#                f[i][k]=[0,0]
#            else:
#                print(f[i][k])
#                f[i][k]=f[i][k]+(d[i][k]-d[j][k])/modulo(d[i][k]-d[j][k])**3*m[i]*m[j]*G
#        a[i][k]=f[i][k]/m[i]

for k in range(nt-1):
    for i in range(n):              # Accelerations
        f[i,k]=[0,0]
        for j in range(n):
            if j==i:
                vaux[j,k]=[0,0]
            else:
                #print(f[i,k])
                vaux[j,k]=(d[j,k]-d[i,k])/modulo(d[j,k]-d[i,k])**3*m[i]*m[j]*G
            #print(i,j,k,vaux[j,k],d[i,k],d[j,k])
            f[i,k]=f[i,k]+vaux[j,k]
            #plotvectorpunto(vaux[j,k]*escala,d[i,k])
        a[i,k]=f[i,k]/m[i]
        v[i,k+1]=v[i,k]+a[i,k]*(t[k+1]-t[k])
        d[i,k+1]=d[i,k]+(v[i,k]+v[i,k+1])/2*(t[k+1]-t[k])+0.5*a[i,k]*(t[k+1]-t[k])**2
        
for i in range(n):
    for j in range(nt):
        x=d[i,j][0]
        y=d[i,j][1]
        plt.plot(x,y,',')
plt.show()

#        print(i,k,f[i,k],d[i,k])
#plt.show()

# for i in range(n):
#     plotvector(d[i,0])
# plt.show()

# for i in range(n):
#     plotvector(f[i,0])
# plt.show()

# for i in range(n):
#     plotvectorpunto(f[i,0]*escala,d[i,0])
#     plotvector(d[i,0])
# plt.show()

#print(t)
#print(v)
#print(a[0][0])
#print(G, ME, MS)

