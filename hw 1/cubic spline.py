# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 14:53:59 2016

@author: Louis
"""
import numpy
import scipy
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

apoints = numpy.transpose([[0,0],[1,1],[2,0],[3,-1],[4,0]])


A =numpy.zeros(((4*(len(apoints[0])-1)),4*(len(apoints[0])-1)))
b =numpy.zeros((4*(len(apoints[0])-1),1))
row = 0

#set up the matrix
for i in range(len(apoints[0])-1):
    #sets the constrain that the polynomial intersects with the lower point    
    for k in range(4): 
        A[row][4*i+k]= apoints[0][i]**k
    b[row] = apoints[1,i]    
    row = row + 1    
    #sets the constrain that the polynomial intersects with the Upper point     
    for k in range(4):
        A[row][4*i+k]= apoints[0][i+1]**k
    b[row] = apoints[1,i+1]
    row = row +1
    #first derivative constraint
    for k in [1,2,3]:
        if(k==1):
           A[row][4*i+k]= apoints[0][i+1]**(k-1)
           if (i<3):
               A[row][4*(i+1)+k] = -apoints[0][i+1]**(k-1)
           else:
               A[row][k] = -apoints[0][0]**(k-1) 
        elif(k==2):
            A[row][4*i+k]= 2*apoints[0][i+1]**(k-1)
            if (i<3):
                A[row][4*(i+1)+k] = -2*apoints[0][i+1]**(k-1)
            else:
                A[row][k] = -2*apoints[0][0]**(k-1)        
        else:
            A[row][4*i+k]= 3*apoints[0][i+1]**(k-1)
            if (i<3):
                A[row][4*(i+1)+k] = -3*apoints[0][i+1]**(k-1)
            else:
                A[row][k] = -3*apoints[0][0]**(k-1)                    
    row = row + 1
    #first derivative constraint
    for k in [2,3]:
        if(k==2):
            A[row][4*i+k]= 2*apoints[0][i+1]**(k-2)
            if (i<3):
                A[row][4*(i+1)+k] = -2*apoints[0][i+1]**(k-2)
            else:
                A[row][k] = -2*apoints[0][0]**(k-2)                 
        else:
            A[row][4*i+k]= 6*apoints[0][i+1]**(k-2)
            if (i<3):
                A[row][4*(i+1)+k] = -6*apoints[0][i+1]**(k-2)
            else:
                A[row][k] = -6*apoints[0][0]**(k-2)                 
    row = row + 1
     
constants = numpy.linalg.solve(A,b)
xvalue = numpy.linspace(0,4,1000)
yvalues= numpy.zeros((1000,1))
for j in range(1000):
    yvalues[j]=(constants[4*(j/250)] +constants[4*(j/250)+1]*(xvalue[j]) + constants[4*(j/250)+2]*(xvalue[j])**2 +constants[4*(j/250)+3]*(xvalue[j])**3)
sxt = scipy.interpolate.CubicSpline(apoints[0],apoints[1],bc_type='periodic')

diff = []
for i in range(1000):
    diff.append(numpy.subtract(yvalues[i],sxt(xvalue[1]*i)))
plt.plot(diff)
#plt.plot(xvalue,yvalues,xvalue,sxt(xvalue))




