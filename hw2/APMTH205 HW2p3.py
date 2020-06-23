"""
Created on Fri Sep 30 16:37:35 2016

@author: Louis
"""
import numpy
import scipy
import matplotlib.pyplot as plt

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
def prob2c(A):
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
                L[i][j] = numpy.abs(U[i][j]/U[j][j])%2
            for k in range(j,n):
                U[i][k] = numpy.abs(U[i][k]-L[i][j]*U[j][k])%2
                
    return(P,L,U)


def GenA(m,n): 
    A = numpy.zeros((m*n , m*n))
    for i in range(n*m):
        A[i][i] = A[i][i]+ 1 # toggles light on push site
        if (i-n>=0):
            A[i][(i-n)] = (A[i][(i-n)]+ 1)%2 #toggles light on sight above push site
        if ((i+n)<(n*m)):
            A[i][(i+n)] = (A[i][(i+n)]+ 1)%2 #toggles light on side below push site
        if((i+1)/n ==i/n):
           A[i][(i+1)] = (A[i][(i+1)]+ 1)%2 # toggles light on side to the right of push site
        if((i-1)/n == i/n):
            A[i][i-1] = (A[i][i-1]+ 1)%2 #toggles light on side to the left of push site
    return(A)
        
            
def nullA():
    report = numpy.zeros((9,9))
    for n in range(1,10):
        for m in range(1,10):
            (P,L,U) = prob2c(GenA(m,n))
            diagonal = 0                
            for i in range(n*m):
                diagonal = diagonal + U[i][i]
            report[m-1][n-1] = n*m-diagonal
            
    return(report)
    
    
    
                  
        
    
    