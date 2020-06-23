# -*- coding: utf-8 -*-
"""
Created on Fri Dec 09 10:09:11 2016

@author: Louis
"""

import numpy
import scipy.sparse
import scipy.integrate
import matplotlib.pyplot as plt
import matplotlib.animation as an
import time
import random

"""This initialized the functions that return the physical constants as a function of temperature (Kelvin).

thermal conductivity in W/(K-m)
specific heat in J/kg
density kg/m^3
dynamic viscosity in Ns/m^2
kinematic viscosity in m^2/s
"""
execfile("physcon.py")


"""initializes the variables
"""

nr=15
nz=20
nt = 1800
tfin = 60 #seconds
zfin = 7.62  #meters
T0=273.15
OD = .375*2.54/100 #m
ID = 5./16 *2.54/100 #m
pipe = numpy.zeros([nr+2,nz+2,nt])
rad = numpy.linspace(0,OD/2,nr+2)
copperbound=0
for i in range(nr)[::-1]:
    if(rad[i] < ID/2.):
        copperbound = i+1
        break
times = numpy.linspace(0,tfin,nt)
distance = numpy.linspace(0,zfin,nz)
dz = distance[1]#meters
dr = rad[1] #meters
dt =times[1] # seconds
M = numpy.zeros([(nr+2)*(nz+2),(nr+2)*(nz+2)])
vm = .54 # m/s 
aw = k(320)[0]/(cp(320)[0]*p(320)[0])
ac = k(320)[2]/(cp(320)[2]*p(320)[2])
picpipe=numpy.zeros([2*nr+1,nz,nt])


def fresh():
    pipe = numpy.zeros([nr+2,nz+2,nt])
    M = numpy.zeros([(nr+2)*(nz+2),(nr+2)*(nz+2)])
   
"""runs the calcualtion for n descrete spatial steps, and nt descrete time steps
"""
def calctube():
    fresh()
    start = time.time()
    for j in range(nr+2):
        if(rad[j]<= ID/2):
            pipe[j,:,0] = 288.7
        else:
            pipe[j,:,0] = numpy.linspace(288.7,373.15,(nr+2-copperbound))[j-copperbound]
  
    tempvec = numpy.reshape(pipe[:,:,0].T,((nz+2)*(nr+2),1))

   
    for i in range(1,nz+1):# i is z
        for j in range(1,copperbound): #j is r
 #Generic matrix entities
            M[j  +(nr+2)*i,j  +(nr+2)*i] = 1.+ 2.*aw*dt/dr/dr                   
            M[j  +(nr+2)*i,j+1  +(nr+2)*i]=-aw*dt/dr*(1./dr + 1./rad[j]/2.)     
            M[j  +(nr+2)*i,j-1  +(nr+2)*i]=-aw*dt/dr*(1./dr - 1./rad[j]/2.)
            M[j  +(nr+2)*i,j+((nr+2))*(i+1)] = -aw*dt*(1./2./dz) + 2*vm*(1-(2*rad[j]/ID)**2)/2*dt/dz
            M[j  +(nr+2)*i,j+((nr+2))*(i-1)] = +aw*dt*(1./2./dz) - 2*vm*(1-(2*rad[j]/ID)**2)/2*dt/dz
            
#Boundary Conditions
   #z=0 conditions j=1:n+1
            M[j ,j] = 1
                 
            
   #z=final j=1:n+1
            M[j  +(nr+2)*(nz+1),j  +(nr+2)*(nz+1)] = 1. + 2.*aw*dt/dr/dr                   
            M[j  +(nr+2)*(nz+1),j+1  +(nr+2)*(nz+1)]=-aw*dt/dr*(1./dr + 1./rad[j]/2.)     
            M[j  +(nr+2)*(nz+1),j-1  +(nr+2)*(nz+1)]=-aw*dt/dr*(1./dr - 1./rad[j]/2.)
            
            #M[j  +(nr+2)*(nz+1),j+((nr+2))*(nz)] = aw*dt*(1./2./dz) - vm/2*dt/dz
             
        for j in range(copperbound,nr+1):
            
#Generic matrix entities
            M[j  +(nr+2)*i,j  +(nr+2)*i] = 1.+ 2.*ac*dt/dr/dr                   
            M[j  +(nr+2)*i,j+1  +(nr+2)*i]=-ac*dt/dr*(1./dr + 1./rad[j]/2.)     
            M[j  +(nr+2)*i,j-1  +(nr+2)*i]=-ac*dt/dr*(1./dr - 1./rad[j]/2.)
            M[j  +(nr+2)*i,j+((nr+2))*(i+1)] = -ac*dt*(1./2./dz)
            M[j  +(nr+2)*i,j+((nr+2))*(i-1)] = +ac*dt*(1./2./dz)
            
#Boundary Conditions
   #z=0 conditions j=1:n+1
            M[j,j] =1               
            
   #z=final j=1:n+1
            M[j  +(nr+2)*(nz+1),j  +(nr+2)*(nz+1)] = 1.+ 2.*ac*dt/dr/dr                   
            M[j  +(nr+2)*(nz+1),j+1  +(nr+2)*(nz+1)]=-ac*dt/dr*(1./dr + 1./rad[j]/2.)     
            M[j  +(nr+2)*(nz+1),j-1  +(nr+2)*(nz+1)]=-ac*dt/dr*(1./dr - 1./rad[j]/2.)
            #M[j  +(nr+2)*(nz+1),j+((nr+2))*((nz))] = 0*ac*dt*(1./2./dz)
#r=0 for  z 1:nz+1
        M[0   +(nr+2)*i,0   +(nr+2)*i]=1+ 2*aw*dt/dr/dr
        M[0   +(nr+2)*i,1   +(nr+2)*i]=-2*aw*dt/dr/dr
        M[(nr+2)*i,((nr+2))*(i+1)] = -aw*dt*(1./2./dz) + 2*vm/2*dt/dz
        M[(nr+2)*i,((nr+2))*(i-1)] = +aw*dt*(1./2./dz) - 2*vm/2*dt/dz

#r=ID condition all but z=0 and z=final
        M[nr+1   +(nr+2)*i,nr+1   +(nr+2)*i]=1
          
            
#z=0

    # r=0 condition
    #M[0,0] =1 + 2*aw/dr/dr               
    #M[0,1]=-2*aw/dr/dr  
    M[0,0]=1
    #r=outer
    M[nr+1,nr+1]=1
#z=final
#r=0
    M[0  +(nr+2)*(nz+1),0   +(nr+2)*(nz+1)]=1+ 2*aw*dt/dr/dr
    M[0  +(nr+2)*(nz+1),1   +(nr+2)*(nz+1)]=-2*aw*dt/dr/dr
    
    
    
#r=OD    
    M[nr+1+(nr+2)*(nz+1),nr+1+(nr+2)*(nz+1)]=1
    
    for n in range(nt-1):
        #print(n)
        timeplus= numpy.linalg.solve(M[:,:],tempvec)
        timeplus[-(nr+2):]=timeplus[-2*(nr+2):-(nr+2)]
        pipe[:,:,n+1] = numpy.reshape(timeplus,((nz+2),(nr+2))).T
        tempvec[:] = timeplus[:]

        
    tottime = time.time()-start
    picpipe[nr,:] = pipe[1,1:nz+1,:]
    for i in range(1,nr+1):
        picpipe[nr-i,:]=pipe[i+1,1:nz+1,:]
        picpipe[nr+i,:]=pipe[i+1,1:nz+1,:]
    
    print(tottime)
    fig = plt.figure()
    a= fig.add_subplot(2,1,1)
    a.set_title('After Interference Effects (60 sec )')
    a=plt.imshow(picpipe[:,:,nt-1]-T0,clim=[0,100],aspect = .15)  
    plt.xticks([0,nz-1],[0,zfin])
    plt.yticks([0,nr],["OD",0])
    
    
    b = fig.add_subplot(2,1,2)
    b.set_title('Before Interference Effects (14.1 sec )')
    b=plt.imshow(picpipe[:,:,423]-T0,clim=[0,100],aspect = .15)
    plt.xticks([0,nz-1],[0,zfin])
    plt.yticks([0,nr],["OD",0])
    plt.colorbar(orientation = 'horizontal') 
    
       
                
            

"""plots the gth radial position. rad[g] will give you the radius in meters
"""        
def plot(g):
    plt.plot(times,pipe[g,:]-273.15)
    plt.xlim(0,5)
    plt.ylim(-5,105)
    plt.xlabel('Time (s)')
    plt.ylabel('Temp (C)')
    plt.legend([round(rad[g],4)*100])

"""plots the average temperature as a function of time
"""

def avgtemp():
    weights = numpy.linspace(0,ID,copperbound)
    avgtemp = numpy.zeros([nt])
    for i in range(nt):
        avgtemp[i] = numpy.sum(numpy.multiply(weights,pipe[0:copperbound,i]))/numpy.sum(weights)
    plt.figure()
    plt.plot(times,avgtemp-T0)
    plt.xlim(0,tfin)
    plt.ylim(0,100)
    plt.ylabel('Temperature (C)')
    plt.xlabel('Time (s)')    