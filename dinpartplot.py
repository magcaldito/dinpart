#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import random

def plotvector(v):
    x=np.array((0,v[0]))
    y=np.array((0,v[1]))
    plt.plot(x,y)
    
def plotvectorpunto(v,punto):
    x=np.array((punto[0],punto[0]+v[0]))
    y=np.array((punto[1],punto[1]+v[1]))
    plt.plot(x,y)

def modulo(v):
    modulo=(v[0]**2+v[1]**2)**0.5
    return modulo


G= 6.674E-11            # Gravity constant 
ME=5.97219E24           # Mass of the earth
MS=1.9891E30            # Mass of the sun
MM=7.3476E22            # Mass of the moon
DES= 1.50E11            # Distance to the sun
DEM=3.844E8             # Distance to the moon
VES= 29784.8            # Velocity earth around sun
VME=1023                # Velocity moon around earth

scmoon=0.00000000001    # Scale f moon
sc=5E-19

n=8                     # no. particles
nt=1000                  # no. time steps
ts=40000                 # time step
t=np.zeros(nt)          # time vector

m=np.zeros(n)           # particle masses vector  


d=np.zeros((n,nt,2))    # position vectors 
v=np.zeros((n,nt,2))    # velocity vectors
a=np.zeros((n,nt,2))    # acceleration vectors
f=np.zeros((n,nt,2))    # force vectors
vaux=np.zeros((n,nt,2))    # force vectors

print('---------------------------------')


for i in range (nt):
    t[i]=i*ts

# for i in range(n):            # Random initial conditions
#     m[i]=random.uniform(8000,12000)*ME
#     d[i][0]=[random.uniform(0.1,0.2)*DES,random.uniform(0.1,0.2)*DES]
#     v[i][0]=[random.uniform(-1,1)*VES,random.uniform(-1,1)*VES]
#     print(d[i][0][1])

for i in range (n):            # Symetric initial conditions
    m[i]=MS/20
    d[i][0]=[DES*np.sin(2*np.pi*i/n),DES*np.cos(2*np.pi*i/n)]
    v[i][0]=[6*VME*(-np.cos(2*np.pi*i/n)),6*VME*np.sin(2*np.pi*i/n)]
    print(d[i,0], v[i,0])


# m[0]=MS                   # Sun, Earth, Moon system initial conditions
# m[1]=ME
# m[2]=MM

# d[0,0]=[0,0]
# d[1,0]=[0,DES]
# d[2,0]=[0,DES+DEM]

# v[0,0]=[0,0]
# v[1,0]=[VES,0]
# v[2,0]=[VES+VME,0]


for k in range(nt-1):
    for i in range(n):              
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

for i in range(n):
    for j in range(nt):
        x=d[i,j][0]
        y=d[i,j][1]
        plt.plot(x,y,',')
        plotvectorpunto(f[i,j]*sc,d[i,j])
plt.show()

