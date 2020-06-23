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

nr=20 # radial steps
nz=40 #longitudinal steps
nt = 300 #time steps
tfin = 100 #seconds end time
zfin = 7.62  #meters end distance
T0=273.15# for C to K
R3 = .5*2.54/100 #m
R2 = .375*2.54/100 #m
R1 = 5./16 *2.54/100 #m
pipe = numpy.zeros([nr+2,nz+2,nt]) #will store temp data
rad = numpy.linspace(0,R3/2,nr+2) # stores radius data
copperbound=0 # index of copper ID
waterbound = 0 # i ndex of Water ID
for i in range(nr)[::-1]:
    if(rad[i] < R1/2.):
        copperbound = i+1
        break
    
for i in range(nr)[::-1]:
    if(rad[i]<R2/2):
        waterbound = i+1
        break
times = numpy.linspace(0,tfin,nt)
distance = numpy.linspace(0,zfin,nz+2)
dz = distance[1]#meters
dr = rad[1] #meters
dt =times[1] # seconds
M = numpy.zeros([(nr+2)*(nz+2),(nr+2)*(nz+2)])
vm = .478*4 # m/s
vw = -.54*.2 #m/s
awater = k(323.15)[0]/(cp(323.15)[0]*p(323.15)[0]) # thermal diffusivity of water
awort =  k(323.15)[1]/(cp(323.15)[1]*p(323.15)[1]) #thermal diffusivity of wort
acopper = k(323.15)[2]/(cp(323.15)[2]*p(323.15)[2]) #thermal diffusivity of copper
#acopper = k(323.15)[3]/(cp(323.15)[3]*p(323.15)[3]) # stainless steel
cpwater = cp(323.15)[0] # specific heat
pwater = p(323.15)[0]#density
cpwort =cp(323.15)[1]#specific heat
pwort = p(323.15)[1] #density
picpipe=numpy.zeros([2*nr+1,nz,nt]) # make it look good by showing the entire pipe not just the radius
heat = numpy.zeros([nt]) # stores heat gain of water vs time
heatcheck = numpy.zeros([nt]) # in the CFC model this looks at heat loss of wort - consistancy check 
vr= numpy.zeros([nr+2,nz+2]) # radial velocity 
vz = numpy.zeros([nr+2,nz+2])#longitudinal velocity
tapwater = 288.7 #K temp of cooling water 288.7C = 60F
KettleVol = 0.0189271 #m^3 = 5 gallons
heat = numpy.zeros([nt])
scale1 =.0000000008 #scale factor for mixing flow (wort)
scale2= .0000000004# scale factor for mixing flow (water)
#scale = 0 # to get back to laminer flow
vortnumberz = 5 #integer number of vortices in z
vortnumberr = 1# integer number of vortices in r
vortwortrr=vortnumberr*2*numpy.pi/(rad[copperbound]-rad[int(copperbound/4)]) #constant a from eq 16 for wort
vortwaterr=vortnumberr*2*numpy.pi/(rad[nr]-rad[waterbound]) # constant a from eq 16 for water
vortz=vortnumberz*2*numpy.pi/(distance[nz+1]-distance[0]) # constant b from equation 16

# creats the velocity map
for i in range(nz+2):
    for j in range(copperbound/4):
        vz[j,i] = vw*(1-rad[j]**2/rad[copperbound-1]**2)
    for j in range(int(copperbound/4),copperbound):
        vz[j,i]=-scale1*vortz*numpy.sin(vortz*distance[i])*numpy.cos(vortwortrr*(rad[j]-rad[int(copperbound/4)]))/rad[j] + vw*(1-rad[j]**2/rad[copperbound-1]**2)
        vr[j,i]=scale1*vortwortrr*numpy.cos(vortz*distance[i])*numpy.sin(vortwortrr*(rad[j]-rad[int(copperbound/4)]))/rad[j]
    for j in range(waterbound,nr+1):
        vz[j,i]=-scale2*vortz*numpy.sin(vortz*distance[i])*numpy.cos(vortwaterr*(rad[j]-rad[waterbound]))/rad[j] + vm*((rad[waterbound]-rad[j])*(rad[j]-rad[nr])/(rad[nr+1]-rad[waterbound])**2)
        vr[j,i]=scale2*vortwaterr*numpy.cos(vortz*distance[i])*numpy.sin(vortwaterr*(rad[j]-rad[waterbound]))/rad[j]

    #vr[j] =vm #constant flow
"""runs the calcualtion for nr, nr descrete spatial steps, and nt descrete time steps
"""
def calctube(outertemp):
    start = time.time()
    #Initial conditions
       
    for j in range(nr+2):
        if(rad[j]<= R1/2):
            pipe[j,:,0] =numpy.linspace(tapwater,outertemp,nz+2)
            pipe[j,-1,0]=outertemp
        elif(rad[j]<=R2/2):
            pipe[j,:,0] = numpy.linspace(outertemp,tapwater,(waterbound-copperbound+2))[j+1-copperbound]
        else:
            pipe[j,:,0]=tapwater
            
    tempvec = numpy.reshape(pipe[:,:,0].T,((nz+2)*(nr+2),1))

    #creates the matrix M to solve the M * T(t+1) = T(t) expression
    for i in range(1,nz+1):# i is z
        for j in range(1,copperbound): #j is r
 #Generic matrix entities
            M[j  +(nr+2)*i,j  +(nr+2)*i] = 1.+ 2.*awort*dt/dr/dr + 2*awort*dt/dz/dz                   
            M[j  +(nr+2)*i,j+1  +(nr+2)*i]=-awort*dt/dr*(1./dr + 1./rad[j]/2.) + vr[j,i]/2*dt/dr
            M[j  +(nr+2)*i,j-1  +(nr+2)*i]=-awort*dt/dr*(1./dr - 1./rad[j]/2.) - vr[j,i]/2*dt/dr
            M[j  +(nr+2)*i,j+((nr+2))*(i+1)] = -awort*dt*(1./dz/dz) + vz[j,i]/2*dt/dz
            M[j  +(nr+2)*i,j+((nr+2))*(i-1)] = -awort*dt*(1./dz/dz) - vz[j,i]/2*dt/dz
            
#Boundary Conditions
   #z=0 conditions j=1:n+1
            M[j  ,j  ] = 1. + 2.*awort*dt/dr/dr                
            M[j  ,j+1 ]=-awort*dt/dr*(1./dr + 1./rad[j]/2.)    
            M[j  ,j-1 ]=-awort*dt/dr*(1./dr - 1./rad[j]/2.)
            
            
   #z=final j=1:n+1
            M[j  +(nr+2)*(nz+1),j  +(nr+2)*(nz+1)] = 1.              
          
        
        for j in range(copperbound,waterbound): #j is r
 #Generic matrix entities
            M[j  +(nr+2)*i,j  +(nr+2)*i] = 1.+ 2.*acopper*dt/dr/dr + 2*acopper*dt/dz/dz                   
            M[j  +(nr+2)*i,j+1  +(nr+2)*i]=-acopper*dt/dr*(1./dr + 1./rad[j]/2.)    
            M[j  +(nr+2)*i,j-1  +(nr+2)*i]=-acopper*dt/dr*(1./dr - 1./rad[j]/2.) 
            M[j  +(nr+2)*i,j+((nr+2))*(i+1)] = -acopper*dt*(1./dz/dz) 
            M[j  +(nr+2)*i,j+((nr+2))*(i-1)] = -acopper*dt*(1./dz/dz) 
            
#Boundary Conditions
   #z=0 conditions j=1:n+1
            M[j ,j] = 1. + 2.*acopper*dt/dr/dr                
            M[j ,j+1]=-acopper*dt/dr*(1./dr + 1./rad[j]/2.)   
            M[j ,j-1]=-acopper*dt/dr*(1./dr - 1./rad[j]/2.)
            
            
   #z=final j=1:n+1
            M[j  +(nr+2)*(nz+1),j  +(nr+2)*(nz+1)] = 1. + 2.*acopper*dt/dr/dr              
            M[j  +(nr+2)*(nz+1),j+1  +(nr+2)*(nz+1)]=-acopper*dt/dr*(1./dr + 1./rad[j]/2.)   
            M[j  +(nr+2)*(nz+1),j-1  +(nr+2)*(nz+1)]=-acopper*dt/dr*(1./dr - 1./rad[j]/2.)
        for j in range(waterbound,nr+1):
            
#Generic matrix entities
            M[j  +(nr+2)*i,j  +(nr+2)*i] = 1.+ 2.*awater*dt/dr/dr +2*awater*dt/dz/dz                  
            M[j  +(nr+2)*i,j+1  +(nr+2)*i]=-awater*dt/dr*(1./dr + 1./rad[j]/2.)  + vr[j,i]/2*dt/dr   
            M[j  +(nr+2)*i,j-1  +(nr+2)*i]=-awater*dt/dr*(1./dr - 1./rad[j]/2.) - vr[j,i]/2*dt/dr
            M[j  +(nr+2)*i,j+((nr+2))*(i+1)] = -awater*dt/dz/dz + vz[j,i]/2*dt/dz
            M[j  +(nr+2)*i,j+((nr+2))*(i-1)] = -awater*dt/dz/dz - vz[j,i]/2*dt/dz
            
#Boundary Conditions
   #z=0 conditions j=1:n+1
            M[j,j] = 1.
                  
                           
            
   #z=final j=1:n+1
            M[j  +(nr+2)*(nz+1),j  +(nr+2)*(nz+1)] = 1.+ 2.*awater*dt/dr/dr              
            M[j  +(nr+2)*(nz+1),j+1  +(nr+2)*(nz+1)]=-awater*dt/dr*(1./dr + 1./rad[j]/2.)      
            M[j  +(nr+2)*(nz+1),j-1  +(nr+2)*(nz+1)]=-awater*dt/dr*(1./dr - 1./rad[j]/2)
#r=0 for  z 1:nz+1
        M[0   +(nr+2)*i,0   +(nr+2)*i]=1+ 2*awort*dt/dr/dr 
        M[0   +(nr+2)*i,1   +(nr+2)*i]=-2*awort*dt/dr/dr    

        

#r=R3 condition all but z=0 and z=final
        M[(nr+1 )   +(nr+2)*i,(nr+1 )   +(nr+2)*i]=1+ 2*awater*dt/dr/dr  
        M[(nr+1 )  +(nr+2)*i,(nr+1 )-1   +(nr+2)*i]=-2*awater*dt/dr/dr 

            
#z=0

    # r=0 condition
    #M[0,0] =1 + 2*aw/dr/dr               
    #M[0,1]=-2*aw/dr/dr  
    M[0,0]=1
    #r=outer
    M[nr+1,nr+1]=1
#z=final
#r=0
    M[0  +(nr+2)*(nz+1),0   +(nr+2)*(nz+1)]=1

    
    
    
#r=OD    
    M[(nr+1)  +(nr+2)*(nz+1),(nr+1)  +(nr+2)*(nz+1)] = 1.+ 2.*awater*dt/dr/dr              
     
    M[(nr+1) +(nr+2)*(nz+1),(nr+1)-1  +(nr+2)*(nz+1)]=-2*awater*dt/dr*(1./dr)
    
    
    #actually runs the simulation
    for n in range(nt-1):
        #print(n)
        timeplus= numpy.linalg.solve(M[:,:],tempvec)
# sink for water        
        timeplus[-(nr+2)+waterbound:]=timeplus[-2*(nr+2)+waterbound:-(nr+2)]
#sink for wort        
        timeplus[0:copperbound]=timeplus[(nr+2):(nr+2)+copperbound]
        pipe[:,:,n+1] = numpy.reshape(timeplus,((nz+2),(nr+2))).T
        # monitor heat
        for j in range(waterbound-1,nr+1):
            heat[n+1] = heat[n+1] + (pipe[j,nz+1,n+1]-288.7)*dt*vz[j,-2]*dr*numpy.pi*2*rad[j]*cpwater*pwater
        for j in range(copperbound):
            heatcheck[n+1] = heat[n+1] + (373.15 - pipe[j,1,n+1])*dt*vz[j,-2]*dr*numpy.pi*2*rad[j]*cpwater*pwater
        
        tempvec[:] = timeplus[:]
        tottime = time.time()-start
        
           
    # fill the pretty pipe matrix    
    picpipe[nr,:] = pipe[1,1:nz+1,:]
    for i in range(1,nr+1):
        picpipe[nr-i,:]=pipe[i+1,1:nz+1,:]
        picpipe[nr+i,:]=pipe[i+1,1:nz+1,:]
        
    print(tottime)
    return(heat)
    

"""plots the gth radial position. rad[g] will give you the radius in meters
"""        
def plot(g,t):
    plt.plot(times,pipe[g,:,t]-273.15)
    plt.xlim(0,5)
    plt.ylim(-5,105)
    plt.xlabel('Time (s)')
    plt.ylabel('Temp (C)')
    plt.legend([round(rad[g],4)*100])

"""Plots the pretty pipe with fluid flow display
"""
def turbplot(t):
    displayvz = numpy.zeros([2*nr+1,nz])
    displayvr = numpy.zeros([2*nr+1,nz])
    for j in range(nr+1):
        displayvz[nr+j,:] =vz[j,1:nz+1]
        displayvr[nr+j,:] =vr[j,1:nz+1]
        displayvr[nr-j,:] =-1*vr[j,1:nz+1]
        displayvz[nr-j,:] =vz[j,1:nz+1]
    plt.figure(figsize = (9,4))
    plt.imshow(picpipe[:,:,t]-273.15,clim = [0,100],aspect = .3)
    plt.quiver(displayvz,displayvr,alpha = .6)
    plt.xticks([0,nz-1],[0,zfin])
    plt.yticks([0,nr],["OD",0])
 