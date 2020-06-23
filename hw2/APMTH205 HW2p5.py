# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 11:20:49 2016

@author: Louis
"""

"""
Created on Fri Sep 30 16:37:35 2016

@author: Louis
"""
import numpy
import scipy
import matplotlib.pyplot as plt
import time

         
def qrfactor(A):
    R =numpy.array(A,dtype = float)
    n= len(A[0])
    m= len(A)
    G = numpy.identity(m)
    Q = numpy.identity(m) #m>n
    for k in range(n):
        for j in range(k+1,m)[::-1]:
            #construct G[j-1][j]
            G = numpy.identity(m)
            a1 = R[j-1][k]
            a2 = R[j][k]
            c=0
            s=0
            if(numpy.abs(a1) > numpy.abs(a2)):
                c=1./numpy.sqrt(1+(a2/a1)**2)
                s=c*a2/a1
            else:
                s = 1./numpy.sqrt(1+(a1/a2)**2)
                c =s*(a1/a2)
                    
            G[j-1][j-1] = c
            G[j][j] = c
            G[j-1][j] = s
            G[j][j-1] = -s
            
            R = numpy.dot(G,R)
            Q = numpy.dot(Q,numpy.transpose(G))
    return(Q,R)
    
def genY(start,stop):
    n=stop-start+1
    Y=numpy.zeros((n,3))
    for i in range(start,stop+1):
            Y[i-start][0]= (i*7./120)**2
            Y[i-start][1]= (i*7./120)
            Y[i-start][2]= 1              
    return(Y)

def genYtot():
    n=20
    Y=numpy.zeros((20,9))
    for i in range(1,21):
        if(i<5):
            Y[i-1][0]= (i*7./120)**2
            Y[i-1][1]= (i*7./120)
            Y[i-1][2]= 1
        elif(i<14):
            Y[i-1][3]= (i*7./120)**2
            Y[i-1][4]= (i*7./120)
            Y[i-1][5]= 1
        else:
            Y[i-1][6]= (i*7./120)**2
            Y[i-1][7]= (i*7./120)
            Y[i-1][8]= 1
    return(Y)
            
    
def timearc(h):
    t=numpy.sqrt(2.*h/9.78)
    h=h*.735
    while(h>0.0000000001):
        t=t+numpy.sqrt(8.*h/9.78)
        h=h*.735    
    return(t)
    
    
        
    
  
    
    
    
    
                  
        
    
    