# -*- coding: utf-8 -*-
"""
Spyder Editor

Created on Wed Oct 12 13:58:45 2016

@author: Louis
"""

import numpy
import scipy.integrate
import matplotlib.pyplot as plt
import matplotlib.animation as an
import time
import random

"""
Problem 2
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
c= 3.43e2 #m/s
h = .366 #m
#dt = h/(2*c)/2.
dt = .0005
w = 100* numpy.pi #1/s
p0 = 10#Pa
times = numpy.linspace(0,1.007,1.007/dt)
speedup = numpy.zeros((200,2))
speedup[3:41,:] =(1,99) 
speedup[41:78,:]=(28,72)
speedup[78:123,:]=(18,76)
speedup[123:160,:]=(28,72)
speedup[160:197,:]=(1,99)


global pierce
global p

def drive(n):
    return(p0*numpy.sin(w*n*dt))


def initialize():
    global p
    p=numpy.zeros((100,200,len(times)))
    p[57:61,15:19,1]=drive(1)

def calculate():
    global p
    initialize()
    
    for i in range(2,int(.3/dt)):
        for k in range(4,197):
            for j in range(int(speedup[k][0]),int(speedup[k][1])):
                p[j,k,i]=2*p[j,k,i-1]-p[j,k,i-2] + c*c*dt*dt/(h*h)*(p[j+1,k,i-1]*(1-pierce[j+1,k]) + p[j,k+1,i-1]*(1-pierce[j,k+1]) +p[j-1,k,i-1]*(1-pierce[j-1,k]) +p[j,k-1,i-1]*(1-pierce[j,k-1])-(4-pierce[j-1,k]-pierce[j+1,k]-pierce[j,k-1]-pierce[j,k+1])*p[j,k,i-1]  )
            
        p[57:61,15:19,i]=drive(i)
        if(i%(len(times)/10)==0):
            print("percent complete:")
            print((i*100)/(len(times)))
    

def video():
    global p
    global anim
    fig=plt.figure()
    ims=[]
    for i in range(len(times)):
        t = str(dt*i)
        im = plt.imshow(p[:,:,i],cmap = 'seismic',clim=(-10,10),animated = True)
        ims.append([im])
    anim = an.ArtistAnimation(fig,ims,interval = 50)
    plt.scatter(building[0],building[1],marker = 's', color = 'black',s = 1 )
    plt.show()
        
        
def plot(i):
    global p
    plt.clf()
    plt.imshow(p[:,:,i],cmap = 'seismic',clim=(-10,10))
    plt.colorbar()
    plt.scatter(building[0],building[1],marker = 's', color = 'black',s = 1 )
    plt.title('time =' + str(dt*i) + 's')
    
    
def buildingpoints():
    global building
    building = []
    for i in range(100):
        for k in range(200):
            if(pierce[i,k]==1):
                building.append([k,i])
    building = numpy.transpose(building)
            
            
def pressureexceeds():
    for i in range(len(times)):
        if(abs(p[35,73,i])>=.001):
            print('C @'+str(dt*i))
            break
    for i in range(len(times)):    
        if(abs(p[61,109,i])>=.001):
            print('G @'+str(dt*i))
            break
    for i in range(len(times)):
        if(abs(p[91,188,i])>=.001):
            print('M @'+str(dt*i))
            break
