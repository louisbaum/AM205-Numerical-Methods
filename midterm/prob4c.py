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



def solve(n):
    y=numpy.zeros((2,n))
    time = numpy.linspace(0,8,n)
    h=time[1]-time[0]
    y[0,0],y[1,0] =fexact(time[0])
    y[0,1],y[1,1] =fexact(time[1])
    for i in range(2,n):
        y[0,i] = y[0,i-1] +3./2*h*(7*numpy.pi*y[1,i-1]*numpy.cos(time[i-1])) - 1./2*h*(7*numpy.pi*y[1,i-2]*numpy.cos(time[i-2]))
        y[1,i] = y[1,i-1] +3./2*h*(-7*numpy.pi*y[0,i-1]*numpy.cos(time[i-1])) - 1./2*h*(-7*numpy.pi*y[0,i-2]*numpy.cos(time[i-2]))
    return(y)
    
    
def fexact(t):
    newu=numpy.sin(7*numpy.pi*numpy.sin(t))
    newv=numpy.cos(7*numpy.pi*numpy.sin(t)) 
    return(newu,newv)
    
    
def errorcomp():
    N = [1000,2000,5000,10000,20000,50000,100000,200000]
    h= numpy.zeros(8)
    error = numpy.zeros(8)
    for i in range(len(N)):
        h[i]=8./(N[i]-1)
        y=solve(N[i])
        u,v = fexact(8)
        error[i] = numpy.linalg.norm([y[0,-1]-u,y[1,-1]-v])
    plt.scatter(h,error)
    ax = plt.gca()
    plt.xlim(-.001,.01)
    plt.ylim(10**(-6),.5)
    ax.set_yscale('log')
    plt.xlabel('h (time)')
    plt.ylabel('absolute error')
    return(error)
    