# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: Louis
"""

import numpy
import scipy.integrate
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.animation as an
import time
import random

"""
Problem 2
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
R=2.
w=5.
L=1.
p=1.
g0 = numpy.zeros([40])
s = numpy.linspace(0,R,251)
g0[20] =1.


def estgrad(b):
    gradest = numpy.zeros([40])
    for i in range(40):
        e= numpy.zeros([40])        
        e[:]=b[:] 
        e[i] = e[i]+1e-7
        fr = r(b)
        fre = r(e)
        gradest[i] = (fre-fr)/(1e-7)
    return(gradest)

    
pi= numpy.pi

def sinchain(s):
    sins = numpy.zeros(20)
    for i in range(1,21):
        sins[i-1] = numpy.sin(i*pi*s/R)
    return(sins)
    
def coschain(s):
    coses = numpy.zeros(20)
    for i in range(1,21):
        coses[i-1] = numpy.cos(i*pi*s/R)
    return(coses)
    

def y(s,b):
    ytot = 0
    sins = sinchain(s)
    ytot = numpy.dot(sins,b[20:])
    return(ytot)

def x(s,b):
    xtot = L*s/R
    sins = sinchain(s)
    xtot= xtot+ numpy.dot(sins,b[:20])
    return(xtot)

def dxds(s,b):
    dxdstot = L/R
    coses = coschain(s)
    coeff = numpy.array(range(1,21))/2.*pi
    dxdstot = dxdstot +numpy.dot(numpy.multiply(coses,coeff),b[:20])
    return(dxdstot)

def dyds(s,b):
    dydstot = 0
    coses = coschain(s)
    coeff = numpy.array(range(1,21))/2.*pi
    dydstot = dydstot +numpy.dot(numpy.multiply(coses,coeff),b[20:])
    return(dydstot)
    

""" returns the action given a vector of parameters b
"""

def r(b):
    mu=2000
    Tintegrand = numpy.zeros(251)
    Vintegrand = numpy.zeros(251)
    Ttot = 0
    Vtot=0
    dydsk =0
    dxdsk= 0
    for i in range(251):
        Tintegrand[i] = p*y(s[i],b)**2*w**2
        dxdsk =dxds(s[i],b)
        dydsk =dyds(s[i],b)
        Vintegrand[i] = mu*(numpy.sqrt(dxdsk**2+dydsk**2)-1)**2
    Ttot = scipy.integrate.simps(Tintegrand,s)
    Vtot = scipy.integrate.simps(Vintegrand,s)
    return(Vtot-Ttot)

""" returns the gradient of the action given a vector of parameters b. grad_i corresponds to the deriviative with respect to the ith component of b
"""

def gradr(b):
    mu=2000
    grad = numpy.zeros([40])
    for k in range(1,21):
        Ctot=0
        Dtot=0
        dy=0
        dx= 0
        Dintegrand = numpy.zeros(251)
        Cintegrand = numpy.zeros(251)
        for i in range(251):
            dx =dxds(s[i],b) 
            dy =dyds(s[i],b)
            Cintegrand[i] = 2*dx*k*pi*numpy.cos(k*pi*s[i]/R)/R*mu*(numpy.sqrt(dx**2+dy**2)-1)/numpy.sqrt(dx**2+dy**2)
            Dintegrand[i] = 2*dy*k*pi*numpy.cos(k*pi*s[i]/R)/R*mu*(numpy.sqrt(dx**2+dy**2)-1)/numpy.sqrt(dx**2+dy**2)  - p*2*w*b[k+19]*w*numpy.sin(k*pi*s[i]/R)**2         
        Ctot = scipy.integrate.simps(Cintegrand,s)
        Dtot = scipy.integrate.simps(Dintegrand,s)
        grad[k-1] = Ctot
        grad[k+19] = Dtot  
    return(grad)
    

    
"""plots the y values of the rope from s=0 to s=2
"""

def ploty(b):
    yvalues = numpy.zeros([251])
    for i in range(251):
        yvalues[i]= y(s[i],b)
    plt.plot(s,yvalues)
    return(yvalues)
            

    

            
            
