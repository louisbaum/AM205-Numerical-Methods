# -*- coding: utf-8 -*-
"""
Created on Thu Nov 03 12:38:45 2016

@author: Louis
"""

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
Problem 3
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
global temp

def calculate():
    global temp
    global times
    global A
    temp[:,0]=1
    temp[52,0] = 1-2*dr/a
    A = numpy.zeros((int(.5/dr+3),int(.5/dr+3)))
    for i in range(1,52):
        r = r1+(i-1)*dr
        A[i,i] = (1.+2.*a/dr)
        A[i,i+1]= -1.*(a/dr+a/2./r)
        A[i,i-1]= -1.*(a/dr-a/2./r)     
        A[0,0]=1
        A[52,52]=1
#    
    # boundry condition enforcement
       
    for i in range(0,len(times)-1):
        temp[:,i+1] = numpy.linalg.solve(A,temp[:,i])
        temp[0,i+1] = temp[2,i+1]
        temp[52,i+1] = temp[50,i+1]- 2*dr/a*temp[51,i+1]
              
        
def plot(j):
    global temp
    global times
    image = numpy.zeros((601,601))
    for i in range(-205,205):
        for k in range(-205,205):
            if(int(numpy.sqrt(i**2+k**2)) > 150 and int(numpy.sqrt(i**2+k**2)) <=200):
                image[i+300,k+300] = temp[int(numpy.sqrt(i**2+k**2)-150),j]
    plt.clf()
    fig = plt.imshow(image, clim=(-1,1),cmap = 'seismic')
    fig.axes.get_xaxis().set_visible(False)
    fig.axes.get_yaxis().set_visible(False)
    plt.title('time =' + str(dt*j) + 's')
    

