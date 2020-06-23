# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 13:58:45 2016

@author: Louis
"""

import numpy
import scipy.integrate
import matplotlib.pyplot as plt
import time
import random

"""
Problem 1
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""

def prob1a():
    
    converg = 8./9.*numpy.pi
    error = numpy.zeros((50,1))    
    for n in range(1,51):# Here I implement the trapazoidal integration for n intervals
        intervals = numpy.linspace(0,numpy.pi/3,n+1)
        fvalues1 = numpy.array(1./(5./4.-numpy.cos(intervals[0:n])))
        fvalues2 = numpy.array(1./(5./4.-numpy.cos(intervals[1:n+1])))
        integral = intervals[1]*numpy.sum(fvalues1+fvalues2)/2.
        error[n-1] = abs(-integral + converg)
    n = range(1,51)
    lim = numpy.zeros((50,1))
    for i in range(1,51): 
        lim[i-1] = (numpy.pi/3/i)**2/36.*1.5396*numpy.pi # 1.5396 is the maximum of (\int_0^{\pi/3}f(x) dx)'' [over 0,\pi/3]
        
    plt.yscale('log')
    plt.xscale('log')
    plt.scatter(range(2,51),error[1::])
    plt.plot(n,lim)
    plt.title("Error vs number of intervals \n Error bound as solid line \n error for n=1 not shown for clarity")
    plt.xlabel("n intervals")
    plt.ylabel("absolute error")
    
def prob1b():
    
    converg = 8./3.*numpy.pi
    error = numpy.zeros((50,1))    
    for n in range(1,51): # trapazoidal integration for n intervals
        intervals = numpy.linspace(0,2*numpy.pi,n+1)
        fvalues1 = numpy.array(1./(5./4.-numpy.cos(intervals[0:n]))) #fk
        fvalues2 = numpy.array(1./(5./4.-numpy.cos(intervals[1:n+1])))#fk+1
        integral = intervals[1]*numpy.sum(fvalues1+fvalues2)/2.
        error[n-1] = abs(-integral + converg)
        
    plt.yscale('log')
    plt.xscale('log')
    plt.scatter(range(1,51),error[0::])
    plt.title("Error vs number of intervals \n error for n=1 not shown for clarity")
    plt.xlabel("n intervals")
    plt.ylabel("absolute error")
    

"""
Problem 2
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Uses Integrate, adaptive and FofX defined below to answer problem 2b
"""
def prob2b():
    global m
    global count
    for m in range(4,9):
        count = 1
        boop =adaptive(-1,.75,10**(-6))
        print(boop)
        print(count)
    


"""
The functions integrate, adaptive and Fofx are coupled. you need to define the function F(x) "Fofx" that you wish to integrate.

integrate integrates F(x) from x=a to x=b
""" 
def integrate(a,b):
    c = (a+b)/2.
    l = (b-a)*1.
    u = numpy.array([-numpy.sqrt(3./5),0,numpy.sqrt(3./5)])
    w = [5./9, 8./9, 5./9]
    x =c+l/2.*u
    f=numpy.zeros((3,1))
    for i in range(3):
        f[i]= Fofx(x[i])
    integrand = l/2.*(numpy.dot(w,f))
    return([integrand[0],f])

"""

Adaptive integrates F(x) from x=a to x=b and if splitting the region into smaller intervals increases accuracy it does so recursivly.
"""

def adaptive(a,b,tol):
    c = (a+b)/2.
    l = (b-a)*1.
    global count
    whole,f = integrate(a,b)
    part1,f1 = integrate(a,c)
    part2,f2 = integrate(c,b)
    even = f1==f2[::-1]
    error = numpy.abs(whole-part1-part2)
    if(not(any(numpy.transpose(even)[0])) and error < tol*(l)):
        return(numpy.array([whole,error]))
    else:
        count = count +1
        return(adaptive(a,c,tol)+adaptive(c,b,tol))
    
   
"""
given an x value returns F(x)
I know globals are not good form but this is "cleaner" then passing m through "prob2b -> adaptive -> integrate-> Fofx"
"""

"""   
Problem 2b F(x) 
def Fofx(x):
    global m
    return(((x*1.)**m-(x*1.)**2+1.))

"""

def Fofx(x):
    return((500*x**6-700*x**4+245*x**2-3)*(numpy.sin(2*numpy.pi*x))**2)
    

def prob2c():
    global count
    count = 1
    boop =adaptive(-1,1,10**(-9))
    print(boop)
    print(count)
    
"""
Problem 5
------------------------------------------------------------------------------------------------------------------------------------------------------------

intersect determines if a line from (x0,y0) to (x1,y1) intersects a circle centered at (px,py) with radius r

As it turns out the perpendicular distance between a circle and a line can be calculated as a dot product of two vectors
1- the vector from an endpoint to the center of the circle.
2 - the unit vector perpendicular to the line segment.

by comparing this perpendicular distance to the radius you can determine if the line intersected the circle.
"""
def intersect(x0,y0,x1,y1,px,py,r):
    vectorA = numpy.array([[(px-x0)*1.],[(py-y0)*1.]]) # vector from point (x0,y0) to the center of the circle
    linevector = numpy.array([[(x1-x0)*1.],[(y1-y0)*1.]])# vector from point (x0,y0) to (x1,y1)
    perpunitvector = numpy.dot([[0.,-1.],[1.,0.]],linevector)/numpy.linalg.norm(linevector)# calculates the unit vector perpendicular to the vector from (x0,y0) to (y0,y1)  
    perpdistance =numpy.abs(numpy.dot(numpy.transpose(vectorA),perpunitvector)[0][0]) # calculates perpendicular distance.
    pardistance = numpy.sqrt(numpy.linalg.norm(vectorA)**2-perpdistance**2)
    return(perpdistance<=r and pardistance <=r)# returns true if it does intersect the circle.
   
    
def asteroid():
    Rearth = 0.02
    Rmoon = 0.005
    Pearth = [0.,0.]
    Pmoon = [1.,0.]
    x0obs = [1.0798,0]
    x1obs = [1.0802,-0.0189]
    n=2500
    vect0 = numpy.zeros((n,4))
    for i in range(n):
        x0 =x0obs[0]+ random.gauss(0,0.002) # random.gauss addes the observation error
        y0 =x0obs[1]+ random.gauss(0,0.002)
        x1 =x1obs[0]+ random.gauss(0,0.002)
        y1 =x1obs[1]+ random.gauss(0,0.002)
        vect0[i,:] =[x0,y0, 1 / 0.02 *(x1-x0),1/0.02*(y1-y0)]
    vect = numpy.zeros((n,501,4))
    t= numpy.linspace(0,10,501)
    for i in range(n):
        vect[i,:,:] = scipy.integrate.odeint(change,vect0[i,:],t)
    
    
    #plots the tragectories
    ax = plt.gca()
    earth = plt.Circle((Pearth[0],Pearth[1]),Rearth,color='green')
    moon = plt.Circle((Pmoon[0],Pmoon[1]),Rmoon, color = 'black')
    for i in range(n):
        plt.plot(vect[i,:,0],vect[i,:,1])
    plt.ylim(-1.2,1.2)
    plt.xlim(-2.4,1)     
    ax.add_artist(earth)
    ax.add_artist(moon)
    plt.show()
    return(vect)
    
        
    
def change(vect,t):# accepts vector(x,y,u,v) returns dvectory/dt
        mu = 0.01
        x,y,u,v = vect
        dvectdt = [u,v, v - x*(1-mu)/(x**2+y**2)**(3./2)+(x-mu)-(x-1)*mu/((x-1)**2+y**2)**(3./2), -u + y - y*(1-mu)/(x**2+y**2)**(3./2)-y*mu/((x-1)**2+y**2)**(3./2) ]
        return(dvectdt)
    
def collide(pos):
    earthind = [] # to keep track of which tragectories collide with the earth
    moonind = []# ditto ---moon
    i=0
    while(i<=2499):
        for k in range(500):
            #check to see if the asteroid collides with the earth
            if(intersect(pos[i,k,0],pos[i,k,1],pos[i,k+1,0],pos[i,k+1,1],0,0,.02)):
                earthind.append([i,k+1])
                break
             #check to see if the asteroid collides with the moon               
            elif(intersect(pos[i,k,0],pos[i,k,1],pos[i,k+1,0],pos[i,k+1,1],1,0,.005)):
                moonind.append([i,k+1])
                break
            
        i=i+1
    return(earthind,moonind)
            
#plots a subset of the tragectories given by the index in list. pos contains x,y,u,v data for all the trajectories.
def samples(pos,list):
    ax = plt.gca()
    earth = plt.Circle((0,0),0.02,color='green')
    moon = plt.Circle((1,0),0.005, color = 'black')
    for i in range(3):
        print(list[i])        
        plt.plot(pos[list[i][0],0:(list[i][1]),0],pos[list[i][0],:(list[i][1]),1])
    plt.ylim(-1.2,1.2)
    plt.xlim(-2.4,1)     
    ax.add_artist(earth)
    ax.add_artist(moon)
    plt.show()
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    