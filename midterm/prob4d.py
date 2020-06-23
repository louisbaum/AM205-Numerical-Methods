# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 11:54:42 2016

@author: Louis
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 10:56:32 2016

@author: Louis
"""
import numpy
import scipy.integrate
import matplotlib.pyplot as plt
import matplotlib.animation as an
import time
import random



def solve(thresh):
    global time
    h=10**(-6)
    ho =[h]
    y=[[],[]]
    time = [0,h]
    u0,v0 =fexact(time[0])
    u1,v1 =fexact(time[1])
    y[0].append(u0)
    y[1].append(v0)
    y[0].append(u1)
    y[1].append(v1)
    t=time[-1]
    i=1
    while(t <= 8):
        test = False
        ht = min(1.1*ho[-1],8-t)
        if(ht==0):
            break
        while(not(test)):
            utest = y[0][i] +(ht*ht+2*ht*ho[-1])/(2*ho[-1])*(7*numpy.pi*y[1][i]*numpy.cos(time[i])) - (ht*ht/(2*ho[-1]))*(7*numpy.pi*y[1][i-1]*numpy.cos(time[i-1]))
            vtest = y[1][i] +(ht*ht+2*ht*ho[-1])/(2*ho[-1])*(-7*numpy.pi*y[0][i]*numpy.cos(time[i])) - (ht*ht/(2*ho[-1]))*(-7*numpy.pi*y[0][i-1]*numpy.cos(time[i-1]))
            yest = ytilde(t,[y[0][i],y[1][i]],ht)
            error = numpy.linalg.norm([yest[0]-utest,yest[1]-vtest])/ht
            test = error <= thresh 
            if(not(test)):
                ht=ht/2   
        ho.append(ht)
        i=i+1        
        time.append(t+ho[-1])
        t = time[-1]
        y[0].append(utest)
        y[1].append(vtest)     
    return(y)
    
    
def fexact(t):
    newu=numpy.sin(7*numpy.pi*numpy.sin(t))
    newv=numpy.cos(7*numpy.pi*numpy.sin(t)) 
    return(newu,newv)
    
    
def errorcomp():
    T = [1.5,2,2.5,3,3.5,4,4.5,5,5.5,6]
    error = numpy.zeros(len(T))
    leny = numpy.zeros(len(T)) 
    for i in range(len(T)):
        y=solve(10**(-1*T[i]))
        u,v = fexact(8)
        error[i] = numpy.linalg.norm([y[0][-1]-u,y[1][-1]-v])
        leny[i] = len(y[0])
        print(T[i], leny[i],error[i])
    return(leny,error)
    
    
def ytilde(t,y,h):
    u=y[0]
    v=y[1]
    K1 = [7*numpy.pi*v*numpy.cos(t),-7*numpy.pi*u*numpy.cos(t)]
    K2 = [7*numpy.pi*(v+K1[1]*h/2)*numpy.cos(t+h/2),-7*numpy.pi*(u+K1[0]*h/2)*numpy.cos(t+h/2)]
    K3 = [7*numpy.pi*(v-K1[1]*h+2*K2[1]*h)*numpy.cos(t+h),-7*numpy.pi*(u-K1[0]*h+2*K2[0]*h)*numpy.cos(t+h)]
    utilde = u+h/6*(K1[0]+4*K2[0]+K3[0]) 
    vtilde = v+h/6*(K1[1]+4*K2[1]+K3[1])
    return([utilde,vtilde])
    
    
    