# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 14:33:44 2016

@author: Louis
"""

import numpy
import scipy
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
n=10000
tvalues = numpy.linspace(0,4,n)
tplusvalues = numpy.linspace(0,4,n)
for i in range(n-1):
    tplusvalues[i] =tplusvalues[i+1]
tplusvalues[n-1] = 0

def ab():
    apoints = numpy.transpose([[0,0],[1,1],[2,0],[3,-1],[4,0]])
    sxt = scipy.interpolate.CubicSpline(apoints[0],apoints[1],bc_type='periodic')
    plt.plot(tvalues, sxt(tvalues),'blue',tvalues,numpy.sin(tvalues * numpy.pi/2),'r')
    plt.title("Sx(t) in blue and cos(tpi/2) in red")
    


def c():
    bpoints = numpy.transpose([[0,1],[1,0],[2,-1],[3,0],[4,1]])
    syt = scipy.interpolate.CubicSpline(bpoints[0],bpoints[1],bc_type='periodic')
    plt.plot(tvalues, syt(tvalues),'green',tvalues,numpy.cos(tvalues * numpy.pi/2),'black')
    plt.title("Sy(t) in green and cos(tpi/2) in black")
    
    
def d(k):
    knum = numpy.linspace(0,2,k)
    kplus = numpy.linspace(0,2,k)
    for i in range(k-1):
        kplus[i] =kplus[i+1]
    kplus[k-1] = 1  
    
    apoints = numpy.transpose([[0,0],[1,1],[2,0],[3,-1],[4,0]])
    sxt = scipy.interpolate.CubicSpline(apoints[0],apoints[1],bc_type='periodic')
    bpoints = numpy.transpose([[0,1],[1,0],[2,-1],[3,0],[4,1]])
    syt = scipy.interpolate.CubicSpline(bpoints[0],bpoints[1],bc_type='periodic')
    plt.scatter(sxt(tvalues),syt(tvalues))

    y=sxt(knum)
    x=syt(knum)
    xup=syt(kplus)
    print(2*numpy.sum(numpy.multiply(y,abs(xup-x))))