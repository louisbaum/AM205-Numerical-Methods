# -*- coding: utf-8 -*-
"""
Created on Sun Dec 11 16:08:15 2016

@author: Louis
"""

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

nr=15# radial points
nz=20 # z points
nt = 18000
tfin = 6000 #seconds
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
distance = numpy.linspace(0,zfin,nz+2)
dz = distance[1]#meters
dr = rad[1] #meters
dt =times[1] # seconds
M = numpy.zeros([(nr+2)*(nz+2),(nr+2)*(nz+2)])

#thermal diffusivity
aw = k(323.15)[0]/(cp(323.15)[0]*p(323.15)[0])
ac = k(323.15)[2]/(cp(323.15)[2]*p(323.15)[2]) # copper
#ac = k(323.15)[3]/(cp(323.15)[3]*p(323.15)[3]) # stainless steel
cpwater = cp(323.15)[0]
pwater = p(323.15)[0]
cpwort =cp(323.15)[1]
pwort = p(323.15)[1]
picpipe=numpy.zeros([2*nr+1,nz,nt])
heat = numpy.zeros([nt])
vr= numpy.zeros([copperbound,nz])
vz = numpy.zeros([copperbound,nz])
tapwater = 285.8 #K
KettleVol = 0.0189271 #m^3 = 5 gallons
heat = numpy.zeros([nt])
scale =.0000003*5
vm = .54*5 # m/s
#scale = 0 # to get back to laminer flow
vortnumberz = 4
vortnumberr = 2
vortr=vortnumberr*2*numpy.pi/(rad[int(copperbound)]-rad[int(copperbound/4)])
vortz=vortnumberz*2*numpy.pi/(distance[nz+1]-distance[1])

# creates the velocity map
for i in range(1,nz+1):
    for j in range(int(copperbound/4),copperbound):
        vz[j,i-1]=-scale*vortz*numpy.sin(vortz*distance[i+1])*numpy.cos(vortr*(rad[j]-rad[int(copperbound/4)]))/rad[j]
        vr[j,i-1]=scale*vortr*numpy.cos(vortz*distance[i+1])*numpy.sin(vortr*(rad[j]-rad[int(copperbound/4)]))/rad[j]
    for j in range(copperbound):
        vz[j,i-1] = vz[j,i-1] + 2.*vm*(1.-rad[j]**2/(ID/2.)**2)
    #vr[j] =vm #constant flow
        
        
"""runs the calcualtion for n descrete spatial steps, and nt descrete time steps
"""
def calctube(outertemp):
    start = time.time()
    for j in range(nr+2):
        if(rad[j]<= ID/2):
            pipe[j,:,0] = tapwater
        else:
            pipe[j,:,0] = numpy.linspace(tapwater,outertemp,(nr+3-copperbound))[j+1-copperbound]

    tempvec = numpy.reshape(pipe[:,:,0].T,((nz+2)*(nr+2),1))

   
    for i in range(1,nz+1):# i is z
        for j in range(1,copperbound): #j is r
 #Generic matrix entities
            M[j  +(nr+2)*i,j  +(nr+2)*i] = 1.+ 2.*aw*dt/dr/dr + 2*aw*dt/dz/dz                   
            M[j  +(nr+2)*i,j+1  +(nr+2)*i]=-aw*dt/dr*(1./dr + 1./rad[j]/2.) + vr[j,i-1]/2*dt/dr    
            M[j  +(nr+2)*i,j-1  +(nr+2)*i]=-aw*dt/dr*(1./dr - 1./rad[j]/2.) - vr[j,i-1]/2*dt/dr
            M[j  +(nr+2)*i,j+((nr+2))*(i+1)] = -aw*dt*(1./dz/dz) + vz[j,i-1]/2*dt/dz
            M[j  +(nr+2)*i,j+((nr+2))*(i-1)] = -aw*dt*(1./dz/dz) - vz[j,i-1]/2*dt/dz
            
#Boundary Conditions
   #z=0 conditions j=1:n+1
            M[j ,j] = 1
            
            
   #z=final j=1:n+1
            M[j  +(nr+2)*(nz+1),j  +(nr+2)*(nz+1)] = 1. + 2.*aw*dt/dr/dr                
            M[j  +(nr+2)*(nz+1),j+1  +(nr+2)*(nz+1)]=-aw*dt/dr*(1./dr + 1./rad[j]/2.)  + vr[j,-1]/2*dt/dr      
            M[j  +(nr+2)*(nz+1),j-1  +(nr+2)*(nz+1)]=-aw*dt/dr*(1./dr - 1./rad[j]/2.)- vr[j,-1]/2*dt/dr
             
        for j in range(copperbound,nr+1):
            
#Generic matrix entities
            M[j  +(nr+2)*i,j  +(nr+2)*i] = 1.+ 2.*ac*dt/dr/dr +2*ac*dt/dz/dz                  
            M[j  +(nr+2)*i,j+1  +(nr+2)*i]=-ac*dt/dr*(1./dr + 1./rad[j]/2.)     
            M[j  +(nr+2)*i,j-1  +(nr+2)*i]=-ac*dt/dr*(1./dr - 1./rad[j]/2.)
            M[j  +(nr+2)*i,j+((nr+2))*(i+1)] = -ac*dt/dz/dz
            M[j  +(nr+2)*i,j+((nr+2))*(i-1)] = -ac*dt/dz/dz
            
#Boundary Conditions
   #z=0 conditions j=1:n+1
            M[j,j] = 1.
                  
                           
            
   #z=final j=1:n+1
            M[j  +(nr+2)*(nz+1),j  +(nr+2)*(nz+1)] = 1.+ 2.*ac*dt/dr/dr              
            M[j  +(nr+2)*(nz+1),j+1  +(nr+2)*(nz+1)]=-ac*dt/dr*(1./dr + 1./rad[j]/2.)      
            M[j  +(nr+2)*(nz+1),j-1  +(nr+2)*(nz+1)]=-ac*dt/dr*(1./dr - 1./rad[j]/2)
#r=0 for  z 1:nz+1
        M[0   +(nr+2)*i,0   +(nr+2)*i]=1+ 2*aw*dt/dr/dr + 2*aw*dt/dz/dz 
        M[0   +(nr+2)*i,1   +(nr+2)*i]=-2*aw*dt/dr/dr 
        M[(nr+2)*i,((nr+2))*(i+1)] = -aw*dt*(1./dz/dz) + vz[0,i-1]/2*dt/dz
        M[(nr+2)*i,((nr+2))*(i-1)] = -aw*dt*(1./dz/dz) - vz[0,i-1]/2*dt/dz
        

#r=ID condition all but z=0 and z=final
        M[nr+1   +(nr+2)*i,nr+1   +(nr+2)*i]=1
        M[j  +(nr+2)*i,j-1  +(nr+2)*i]=-ac*dt/dr*(1./dr - 1./rad[j]/2.) 
            
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
         
        for j in range(0,copperbound):
            heat[n+1] = heat[n+1] + (pipe[j,nz+1,n+1]-tapwater)*dt*vz[j,-2]*dr*numpy.pi*2*rad[j]*cpwater*pwater
        for i in range(nz+1):
                timeplus[nr+1+(nr+2)*i] = timeplus[nr+1+(nr+2)*i] - heat[n+1]/(cpwort*pwort*KettleVol)    
        tempvec[:] = timeplus[:]
        tottime = time.time()-start
        
           
        
    picpipe[nr,:] = pipe[1,1:nz+1,:]
    for i in range(1,nr+1):
        picpipe[nr-i,:]=pipe[i+1,1:nz+1,:]
        picpipe[nr+i,:]=pipe[i+1,1:nz+1,:]
        
    print(tottime)
    return(heat)              
            
"""plots the pipe and fluid flow at time t
"""
def turbplot(t):
    displayvz = numpy.zeros([2*nr-1,nz])
    displayvr = numpy.zeros([2*nr-1,nz])
    for j in range(copperbound):
        displayvz[nr+j,:] =vz[j,:]
        displayvr[nr+j,:] =vr[j,:]
        displayvr[nr-j,:] =-1*vr[j,:]
        displayvz[nr-j,:] =vz[j,:]
    plt.figure(figsize = (9,4))
    plt.imshow(picpipe[1:2*nr,:,t]-273.15,clim = [0,100],aspect = .3)
    plt.quiver(displayvz,displayvr)
    plt.xticks([0,nz-1],[0,zfin])
    plt.yticks([0,nr],["OD",0])

