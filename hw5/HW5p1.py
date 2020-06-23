# -*- coding: utf-8 -*-
"""
Spyder Editor

Created on Wed Oct 12 13:58:45 2016

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
Problem 1
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""

def function(x0):
    return(100*(x0[1]-x0[0]**2)**2+(1-x0[0])**2)

def grad(x0):
    y=x0[1]
    x=x0[0]
    dx = 2*(200*x**3-200*x*y+x-1)
    dy =200*(y-x**2)
    #boop = numpy.linalg.norm([dx,dy])
    return([dx,dy])

global a
global b
x = numpy.linspace(-2,2,500)
y = numpy.linspace(-2,2,500)
z = numpy.zeros((500,500))
points = [[[-1,1]],[[0,1]],[[2,1]]]
its = [0,0,0]
path = [[[-1,1]],[[0,1]],[[2,1]]]  
step =1

def zvalues():
    for i in range(len(x)):
        for k in range(len(y)):
            z[i,k] = function(numpy.array([x[i],y[k]]))


""" this is part of my line search
"""
def minfunc(x0):
    x=x0[0]
    y=x0[1]
    dx,dy = grad([x,y])
    x = x0[0] - x0[2]*dx
    y = x0[1] - x0[2]*dy
    return(function([x,y]))
    

"""
Here is my script for steepest descent
"""
def steepdescent():
    for i in range(len(points)):        
        step = 1
        x=points[i][0][0]
        y=points[i][0][1]        
        while(its[i] <= 2000 and step >= 1e-8 ):        
            dx,dy = grad([x,y])
            guess = 1/(dx**2+dy**2)**(.5)/(its[i]+1)# here modify the initial guess based on iteration - the further along we are the shorter the step we expect
            eta = scipy.optimize.minimize(minfunc,[x,y,guess],bounds=((x,x),(y,y),(0,None))) # here we implement a lineserch
            y = y-eta.x[2]*dy
            x = x-eta.x[2]*dx
            step = eta.x[2]*numpy.linalg.norm([dx,dy])
            path[i].append([x,y])
            its[i]=its[i]+1
            
"""
Here is my script for newtons method
"""
            
def NM():
    for i in range(len(points)):
        step = 1
        x=points[i][0][0]
        y=points[i][0][1]
        while(step >= 1e-8 ):
            H = Hessian(x,y)
            gr = numpy.multiply(-1,grad([x,y]))
            s = numpy.linalg.solve(H,gr)
            x=x+s[0]
            y=y+s[1]
            path[i].append([x,y])
            step = numpy.linalg.norm(s)
            its[i]=its[i]+1
            
"""
Here is my script for BFGS
"""
    
            
def BFGS():
    for i in range(len(points)):
        step = 1
        H=numpy.identity(2)
        x=points[i][0][0]
        y=points[i][0][1]
        while(step >= 1e-8 ):
            gr = numpy.multiply(-1,grad([x,y]))
            s = numpy.dot(H,gr)
            x=x+s[0]
            y = y+s[1]
            grup = grad([x,y])
            boop = numpy.add(grup,gr)
            path[i].append([x,y])
            step = numpy.linalg.norm(s)
            p = 1./(numpy.dot(boop,s))
            A = numpy.add(numpy.identity(2),[[-1*p*s[0]*boop[0],-1*p*s[0]*boop[1]],[-1*p*s[1]*boop[0],-1*p*s[1]*boop[1]]])
            B =numpy.add(numpy.identity(2),[[-1*p*s[0]*boop[0],-1*p*s[1]*boop[0]],[-1*p*s[0]*boop[1],-1*p*s[1]*boop[1]]])
            C = [[p*s[0]*s[0],p*s[0]*s[1]],[p*s[1]*s[0],p*s[1]*s[1]]]
            H = numpy.add(numpy.dot(numpy.dot(A,H),B),C)
            its[i]=its[i]+1       



def Hessian(x,y):
    H00=1200*x**2-400*y+2
    H01=-400*x
    H10=-400*x
    H11=200
    return([[H00,H01],[H10,H11]])


"""These function plot the paths from the three points
"""        
def plot1(title):
    plt.contour(x,y,z,[0.1,1,10,50,100,250,500,1000],alpha = .6)
    plt.plot(numpy.transpose(path[0])[0],numpy.transpose(path[0])[1],linewidth = 2, color = 'black')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(title+'\n''Path from [-1,1]')

def plot2(title):
    plt.contour(x,y,z,[0.1,1,10,50,100,250,500,1000],alpha = .6)
    plt.plot(numpy.transpose(path[1])[0],numpy.transpose(path[1])[1],linewidth = 2, color = 'black')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(title+'\n''Path from [0,1]')
    
def plot3(title):
    plt.contour(x,y,z,[0.1,1,10,50,100,250,500,1000],alpha = .6)
    plt.plot(numpy.transpose(path[2])[0],numpy.transpose(path[2])[1],linewidth = 2, color = 'black')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(title+'\n''Path from [2,1]')
    
            
            
