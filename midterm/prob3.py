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

dx = 10**(-5)

def svd(F):
    u,s,v = numpy.linalg.svd(F)
    sigma  =[[s[0],0],
             [0,s[1]]]
    R = numpy.dot(u,v)
    P = numpy.dot(numpy.transpose(v),numpy.dot(sigma,v))
    return(R,P)   
    
def genF(Xvec):
    X=Xvec[0]
    Y=Xvec[1]
    F =[[1+.2*numpy.cos(3*Y)*numpy.cos(X)-.1*numpy.sin(X+Y),-.3*numpy.sin(X)*numpy.sin(3*Y)-0.1*numpy.sin(X+Y)],[-.3*numpy.sin(3*X+2*Y),1-0.2*numpy.sin(3*X+2*Y)]]
    return(F)


"""
determines F, R, P in teh neighborhood of the point Xvec and returns the gradient calculed with the 4th order accurate equation found in prob 2
"""
    
def gradG(Xvec):
    X=Xvec[0]
    Y=Xvec[1]
    xpos = [X-3*dx,X-2*dx,0,X,0,X+2*dx,X+3*dx]
    ypos = [Y-3*dx,Y-2*dx,0,Y,0,Y+2*dx,Y+3*dx]
    Gx = numpy.zeros((7,3))
    Gy = numpy.zeros((7,3))
    gradG = numpy.zeros((2,3))
    for i in range(7):
        F=genF([xpos[i],Y])
        R,P=svd(F)
        Gx[i,0]=R[0,1]
        Gx[i,1]=numpy.linalg.det(F)
        Gx[i,2]=numpy.linalg.norm(numpy.subtract(P,numpy.identity(2)))  
    for k in range(7):
        F=genF([X,ypos[k]])
        R,P=svd(F)
        Gy[k,0]=R[0,1]
        Gy[k,1]=numpy.linalg.det(F)
        Gy[k,2]=numpy.linalg.norm(numpy.subtract(P,numpy.identity(2)))
    for j in range(3):
        gradG[0,j] = (2./15*Gx[0,j]-9./20*Gx[1,j]+9./20*Gx[5,j]-2./15*Gx[6,j])*dx
        gradG[1,j] = (2./15*Gy[0,j]-9./20*Gy[1,j]+9./20*Gy[5,j]-2./15*Gy[6,j])*dx
    return(gradG)

def G(Xvec):
    G= numpy.zeros(3)
    F=genF([Xvec[0],Xvec[1]])
    R,P=svd(F)
    G[0]=R[0,1]
    G[1]=numpy.linalg.det(F)
    G[2]=numpy.linalg.norm(numpy.subtract(P,numpy.identity(2)))
    return(G)
                     
def findzero(Xvec,i):
    X0=Xvec[0]
    Y0=Xvec[1]
    hessian = numpy.zeros((2,2))
    error=1
    dx = .5*10**(-6)   
    while( error > 10**(-6)):
        grG= gradG([X0,Y0])
        grGxplus = gradG([X0+dx,Y0])
        grGxminus = gradG([X0-dx,Y0])
        grGyplus = gradG([X0,Y0+dx])
        grGyminus = gradG([X0,Y0-dx])
        delf = numpy.transpose([grG[0,i],grG[1,i]])
        hessian[0,:] =[(grGxplus[0,i]-grGxminus[0,i])/(2*dx),(grGyplus[0,i]-grGyminus[0,i])/(2*dx)]
        hessian[1,:]= [(grGxplus[1,i]-grGxminus[1,i])/(2*dx),(grGyplus[1,i]-grGyminus[1,i])/(2*dx)]
        # This if chain adapts the step size to the gradient so that this process doesnt take forever.
        dX = -0.1*numpy.dot(numpy.linalg.inv(hessian),delf)[0]
        dY = -0.1*numpy.dot(numpy.linalg.inv(hessian),delf)[1]
        X1 = X0  +dX
        Y1 = Y0 + dY
        error = numpy.linalg.norm([dX, dY])*100
        X0=X1
        Y0=Y1
        if(abs(Y1)>numpy.pi or abs(X1)>numpy.pi):
            return([0,0])
    Xf=X1
    Yf=Y1    
    return([Xf,Yf])
    
def findallzeros(n,i):
    grid = numpy.linspace(-numpy.pi,numpy.pi,n)
    dup = []
    for j in grid:
        for k in grid:
            candidate = findzero([j,k],i)
            dup.append(candidate)
    extrema = []
    return(dup)        
    
def unique(list):
   dummy = [a for a in list if a != [0.,0.] ]
   dummy = numpy.multiply(dummy,10**5)
   dummy= numpy.trunc(dummy)
   print(len(dummy))
   uni = numpy.zeros((len(dummy),2))
   
   for i in range(len(dummy)):
       flag = False
       for k in range(len(uni)):
           if(all(dummy[i,:]==uni[k])):
               flag = True
       if(flag == False):
           uni[i,:] = dummy[i,:]
   uni =numpy.multiply(numpy.asarray([a for a in uni if all(a != [0.,0.])]),10**(-5))
   print(len(uni))
   return(uni)
   
def crittype(pts,i):
    hessian = numpy.zeros((2,2))
    crittype = []
    for j in pts:
        X0=j[0]
        Y0=j[1]
        grG= gradG([X0,Y0])
        grGxplus = gradG([X0+dx,Y0])
        grGxminus = gradG([X0-dx,Y0])
        grGyplus = gradG([X0,Y0+dx])
        grGyminus = gradG([X0,Y0-dx])
        delf = numpy.transpose([grG[0,i],grG[1,i]])
        hessian[0,:] =[(grGxplus[0,i]-grGxminus[0,i])/(2*dx),(grGyplus[0,i]-grGyminus[0,i])/(2*dx)]
        hessian[1,:]= [(grGxplus[1,i]-grGxminus[1,i])/(2*dx),(grGyplus[1,i]-grGyminus[1,i])/(2*dx)] 
        eig = numpy.linalg.eigvals(hessian)
        pos=0
        neg=0
        for k in eig:
            if (k>0):
                pos= pos+1
            elif (k<0):
                neg = neg+1
        if (neg ==0 and pos==2):
            crittype.append(1)
        elif (neg ==2 and pos ==0):
            crittype.append(-1)
        elif(neg==1 and pos ==1):
            crittype.append(0)
    return(crittype)
    
def findmax(list):
    Gval = numpy.zeros((3,len(list)))
    k=0
    for i in list:
        dum = numpy.transpose(G(i))
        Gval[0,k]= dum[0]
        Gval[1,k]=dum[1]
        Gval[2,k]=dum[2]
        k=k+1
    g1max= numpy.argmax(Gval[0])   
    g2max = numpy.argmax(Gval[1]) 
    g3max = numpy.argmax(Gval[2])
    print(list[g1max])
    print(list[g2max])
    print(list[g3max])

def findmin(list):
    Gval = numpy.zeros((3,len(list)))
    k=0
    for i in list:
        dum = numpy.transpose(G(i))
        Gval[0,k]= dum[0]
        Gval[1,k]=dum[1]
        Gval[2,k]=dum[2]
        k=k+1
    g1max= numpy.argmin(Gval[0])   
    g2max = numpy.argmin(Gval[1]) 
    g3max = numpy.argmin(Gval[2])
    print(list[g1max])
    print(list[g2max])
    print(list[g3max])      
    
def showpoints(pts,crittype):
        points = numpy.transpose(pts)
        plt.scatter(points[0],points[1],c=crittype,s = 50,clim= [-2,2],cmap = 'bwr')