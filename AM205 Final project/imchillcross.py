# -*- coding: utf-8 -*-
"""
Created on Fri Dec 09 10:04:31 2016

@author: Louis
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Dec 05 14:54:02 2016

@author: Louis 
"""

import numpy
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
n=10 # number of radial steps
nt = 600 # number of time steps
tfin = 60 # end time
T0=273.15 # for K to C transitions
OD = .375*2.54/100 #R2 m
ID = 5./16 *2.54/100 # R1m
pipe = numpy.zeros([n+2,nt])
M = numpy.zeros([n+2,n+2])
rad = numpy.linspace(0,OD/2,n+2)
copperbound=0
for i in range(n)[::-1]:
    if(rad[i] < ID/2.):
        copperbound = i+1
        break
times = numpy.linspace(0,tfin,nt)
dr = rad[1]-rad[0]
dt =times[1] # seconds

def fresh():
    pipe = numpy.zeros([n+2,nt])
    M = numpy.zeros([n+2,n+2])
   

"""runs the calcualtion for n descrete spatial steps, and nt descrete time steps
"""
def calc():
    fresh()
    for i in range(n+2):
        if(rad[i]<= ID/2):
            pipe[i,0] = 288.7
        else:
            pipe[i,0] = numpy.linspace(288.7,373.15,(n+2-copperbound))[i-copperbound]
            
    for j in range(nt-1):
        if(j%(nt/10)==0):
            print(str((10*j)/nt) + "0% complete")
        for i in range(1,n+1):
            if(rad[i]<= ID/2):
                tk = k(pipe[i,j])[0]# 0 for water
                tcp= cp(pipe[i,j])[0]
                tp=  p(pipe[i,j])[0]
            else:
                tk = k(pipe[i,j])[2] # 2 for copper
                tcp= cp(pipe[i,j])[2]
                tp=  p(pipe[i,j])[2]
            a = 1.*tk/tcp/tp*dt
            M[i,i] = 1. + 2.*a/dr/dr
            M[i,i+1]=-a/dr*(1./dr + 1./rad[i]/2.)
            M[i,i-1]=-a/dr*(1./dr - 1./rad[i]/2.)
        M[0,0]=1+ 2*a*dt/dr/dr
        M[0,1]=-2*a*dt/dr/dr
        M[n+1,n+1]=1
        pipe[:,j+1] =  numpy.linalg.solve(M,pipe[:,j])

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