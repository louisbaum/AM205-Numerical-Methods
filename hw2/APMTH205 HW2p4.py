"""
Created on Fri Sep 30 16:37:35 2016

@author: Louis
"""
import numpy
import scipy
import matplotlib.pyplot as plt
import time

"""
Part A)
 to check for a 0 diagonal element it takes the product of the diagonals and determines it the result is 0
"""

def fsolve(L,b):
    check = 1
    n = len(L)
    for i in range(n):
        check = check*L[i][i] #product of diagonal elements
    if(check == 1):
        x = numpy.zeros((n,1))
        for i in range(n):  #calculates x
            x[i] =(b[i]-numpy.dot(numpy.transpose(x),L[i]))/ L[i][i]       
        return(numpy.mod(x,2))  
        
    else:
        print("L is a Singular Matrix")
        return(None)
 
"""
Part B)
 to check for a 0 diagonal element it takes the product of the diagonals and determines it the result is 0
"""       
def rsolve(U,b):
    check = 1
    n = len(U)
    for i in range(n):
        check = check*U[i][i] #product of diagonal elements
    if(check == 1):
        x = numpy.zeros((n,1))
        for i in range(n)[::-1]: #calculates x
            x[i] =(b[i]-numpy.dot(numpy.transpose(x),U[i]))/ U[i][i]      
        return(numpy.mod(x,2))  
        
    else:
        print("U is a Singular Matrix")
        return(None)
    
"""
This literally implements the pseudocode from the lecture 7 slide


I have ammended this to work only on binary matricies
"""    
def prob2cnotb(A):
    n = len(A)
    P = numpy.identity(n)
    L = numpy.identity(n)
    U = numpy.array(A,dtype = float)
    for j in range(n):
        pivot =  numpy.argmax(numpy.transpose(U)[j][j:n])+j
        U[[j,pivot],j:n] = U[[pivot,j],j:n]
        L[[j,pivot],0:j] = L[[pivot,j],0:j]
        P[[j,pivot],:] = P[[pivot,j],:]
        for i in range(j+1,n):
            if(U[j][j] ==0):
                L[i][j]=0
            else:
                L[i][j] = U[i][j]/U[j][j]
            for k in range(j,n):
                U[i][k] = U[i][k]-L[i][j]*U[j][k]
                
    return(P,L,U)


    
"""
Generates a matrix with ones of the diagonal and -1 below the diagonal.

adds a lower triagonal matrix of -1 to 2* the identity matrix
"""
def generate_g(n):
    Matrix = numpy.zeros((n,n))
    for i in range(n):
        Matrix[i][i] =1
        Matrix[i][n-1] =1
        for k in range(i+1,n):
            Matrix[k][i] = -1
            
    return(Matrix)
    
    
def timegen():
    n = range(10,1000,10)
    t = numpy.zeros((99,1))
    for i in range(99):
        tstart = time.time()
        generate_g(n[i])
        t[i] = time.time()-tstart
    plt.scatter(n,t)
    return(n,t)
    
def residual():
    n = range(10,1000,10)
    error = numpy.zeros((99,1))
    for i in range(99):
        G = generate_g(n[i])
        x = numpy.ones((n[i],1))
        b= numpy.dot(G,x)
        xnew = numpy.linalg.solve(G,b)
        error[i] = numpy.linalg.norm(numpy.subtract(x,xnew))*1./numpy.linalg.norm(xnew)
    plt.scatter(n,error)
        
    
    
  
    
    
    
    
                  
        
    
    