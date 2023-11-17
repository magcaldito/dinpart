#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import random
import pygame

pygame.init()
# pygame.font.init()
# pygame.mixer.init()

white=pygame.Color(255,255,255)
green=pygame.Color(0,255,0)        
yellow=pygame.Color(255,255,0)      
cyan=pygame.Color(0,255,255)      
blue=pygame.Color(0,0,255)
black=pygame.Color(0,0,0)
pink=pygame.Color(255,100,100)
grey=pygame.Color(50,50,50)
ANCHOP=1020
ALTOP=710

pygame.display.set_caption('Particle Dynamics')
ventana= pygame.display.set_mode((ANCHOP,ALTOP))

G= 6.674E-11            # Gravity constant 
ME=5.97219E24           # Mass of the earth
MS=1.9891E30            # Mass of the sun
MM=7.3476E22            # Mass of the moon
DES= 1.50E11            # Distance to the sun
DEM=3.844E8             # Distance to the moon
VES= 29784.8            # Velocity earth around sun
VME=1023                # Velocity moon around earth
DNS=1410/1000           # Sun density
DNE=5510                # Earth density
DNM=3350                # Moon density

xmax=1.5*DES                # x coordinate left screen
xc=-0                 # x coordinate center screen  
yc= 0                 # y coordinate center screen

centro=np.array([xc,yc])
centrop=np.array([ANCHOP/2,ALTOP/2])
escala=ANCHOP/2/(xmax-xc)
escala_r=1              # radius scale
escala_v=1E6            # velocity vector scale
escala_f=5E-13          # force vector scale
escala_a=1E12           # acceleration vecto scale
salto=3                 # gap between displayed steps

def dibujavector(v,color):
    x=np.array((0,v[0]))
    y=np.array((0,v[1]))        
    pygame.draw.line(ventana,color,(0,v[0],0,v[1]),width=1)

def dibujavectorpunto(v,punto,color):
    dpunto1= centrop+(punto+centro)*escala
    dpunto2= centrop+(punto+v+centro)*escala
    pygame.draw.line(ventana,color,dpunto1,dpunto2,width=1)

def dibujacirculo(punto,radio,color):
    #dpunto=np.zeros(2)
    dpunto=centrop+(punto+centro)*escala
    dradio=int(radio*escala)
    if dradio < 1:
        dradio =1
    pygame.draw.circle(ventana,color,dpunto,dradio,width=1)

def modulo(v):
    modulo=(v[0]**2+v[1]**2)**0.5
    return modulo

def radio(masa,densidad):
    radio=((masa/densidad*3/4/np.pi)**(1/3))
    return radio

n=3                     # no. particles
nt=2000                 # no. time steps
ts=15000                # time step
t=np.zeros(nt)          # time vector

m=np.zeros(n)           # particle masses array  
dn=np.zeros(n)          # particle densities array

d=np.zeros((n,nt,2))    # position vectors 
v=np.zeros((n,nt,2))    # velocity vectors
a=np.zeros((n,nt,2))    # acceleration vectors
f=np.zeros((n,nt,2))    # force vectors
vaux=np.zeros((n,nt,2)) # force vectors

print('---------------------------------')


for i in range (nt):
    t[i]=i*ts

# for i in range(n):          # Condciones iniciales aleatorias
#     m[i]=random.uniform(8000,12000)*ME
#     d[i][0]=[random.uniform(0.1,0.2)*DES,random.uniform(0.1,0.2)*DES]
#     v[i][0]=[random.uniform(-1,1)*VES,random.uniform(-1,1)*VES]
#     print(d[i][0][1])


# for i in range (n):            # Symetric initial conditions
#     m[i]=MS/20
#     d[i][0]=[DES*np.sin(2*np.pi*i/n),DES*np.cos(2*np.pi*i/n)]
#     v[i][0]=[6*VME*(-np.cos(2*np.pi*i/n)),6*VME*np.sin(2*np.pi*i/n)]
#     dn[i]=DNE
#     #print(d[i,0], v[i,0])


m[0]=MS                   # Condiciones iniciales Sistema Sol, Tierra, Luna
m[1]=ME
m[2]=MM

d[0,0]=[0,0]
d[1,0]=[0,DES]
d[2,0]=[0,DES+DEM]

v[0,0]=[0,0]
v[1,0]=[VES,0]
v[2,0]=[VES+VME,0]

dn[0]=DNS
dn[1]=DNE
dn[2]=DNM


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
        if k % salto == 1 and k>salto:
            dibujacirculo(d[i,k],radio(m[i],dn[i])*escala_r,pink)
            dibujavectorpunto(f[i,k]*escala_f,d[i,k],yellow)
            dibujavectorpunto(v[i,k]*escala_v,d[i,k],cyan)
            dibujacirculo(d[i,k-salto],radio(m[i],dn[i])*escala_r,black)
            dibujavectorpunto(f[i,k-salto]*escala_f,d[i,k-salto],black)
            dibujavectorpunto(v[i,k-salto]*escala_v,d[i,k-salto],black)
            dibujacirculo(d[i,k-salto],radio(m[i],dn[i])*escala_r/10,grey)
            
        pygame.display.flip()
        
print(radio(ME,DNE))
print(radio(MS,DNS))