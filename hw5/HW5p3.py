# -*- coding: utf-8 -*-
"""
Created on Thu Dec 01 09:56:42 2016

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
Problem 3 - test to see if it works
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""

def vtest(x):
    return(x**2/10.)

def v1(x):
    return(numpy.abs(x))

def v2(x):
    return(12.*(x/10.)**4.-x**2./18. +x/8.+13./10.)  

def v3(x):
    return(8.*numpy.abs(numpy.abs(numpy.abs(x)-1)-1))

"""Modify the function v(x) to return the appropriate potential function.
"""
    
def v(x):
    return(vtest(x))

xpoints = numpy.linspace(-12,12,1921)
h = xpoints[1]-xpoints[0]    
M = numpy.zeros([1923,1923])


"""calculates the appropriate Matrix  and the eigen system for that matrix.
"""
def calc():
    for i in range(1,1922):
        M[i,i] = 2./(h)**2 + v(xpoints[i-1])
        M[i,i+1] = -1/(h)**2
        M[i,i-1] = -1/(h)**2
    
    M[0,0] = 2/h**2 +v(-12)
    M[0,1] = 0

    M[1922,1922] = 2/h**2 +v(12)
    M[1922,1921] = 0

    eigs = numpy.linalg.eig(M)
    E=numpy.sort(eigs[0])[:5]
    indexE = numpy.argsort(eigs[0])[:5]
    return(eigs,indexE,E)

""" Plots the first 5 eigenvalues s 3*Psi_k + E_k
"""

def display(eigs,indexE,E):
    plt.figure()
    plt.plot(xpoints,v(xpoints))
    for i in range(5):
        plt.plot(xpoints,3*eigs[1][1:1922,indexE[i]]+eigs[0][indexE[i]])
    for i in range(5):
        plt.plot([-12,12],[eigs[0][indexE[i]]]*2,color = 'black', alpha = .6)
    #label = ['v(x)','E1','E2','E3','E4','E5']
    #plt.legend(label)
    plt.xlim(-12,12)
    plt.ylim(0,1+numpy.ndarray.max(5*eigs[1][1:1922,indexE[4]]+eigs[0][indexE[4]]))


"""returns the probabilities of the first 5 wavefunctions to be found in the interval [0,6]
"""
    
def prob(param):
    eigs= param[0]
    indexE = param[1]
    E = param[2]
    probs=numpy.zeros([5])
    for i in range(len(E)):
        psivals = numpy.zeros([1921])
        norm=0
        for k in range(len(xpoints)):
                psivals[k]=eigs[1][1+k,indexE[i]]**2
        norm = scipy.integrate.simps(psivals,xpoints)
        probEi = scipy.integrate.simps(psivals[960:1441],xpoints[960:1441])
        probs[i] = probEi/norm
    return(E,probs)

