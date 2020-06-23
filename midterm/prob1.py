# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 09:30:33 2016

@author: Louis
"""

import numpy
import scipy.integrate
import matplotlib.pyplot as plt
import matplotlib.animation as an
import time
import random

F =[[16.,8.,2.,1.],
     [0.,16.,4.,2.],
     [0.,0.,16.,8.],
     [1.,0.,0.,16.]]

end= 20     
Fk=numpy.zeros((4,4,end))
Fk[:,:,0] = F

def iterate():
    for i in range(1,end):
        
        B = numpy.transpose(numpy.linalg.inv(Fk[:,:,i-1]))
        Fk[:,:,i] = .5*numpy.add(Fk[:,:,i-1],B)
    
    print('R = ')
    R = Fk[:,:,end-1]
    print(R)
    
    print('P = ')
    P=numpy.dot(numpy.transpose(Fk[:,:,end-1]),F)
    print(numpy.dot(numpy.transpose(Fk[:,:,end-1]),F))
    
    print('F = ')
    print(numpy.dot(Fk[:,:,end-1],P))
    return(R,P)
    

def svd():
    u,s,v = numpy.linalg.svd(F)
    sigma  =[[s[0],0,0,0],
             [0,s[1],0,0],
             [0,0,s[2],0],
             [0,0,0,s[3]]]
    print('F = ')
    print(numpy.dot(u,numpy.dot(sigma,v)))
    print('R = ')
    R = numpy.dot(u,v)
    print(R)
    
    print('P = ')
    P = numpy.dot(numpy.transpose(v),numpy.dot(sigma,v))
    print(P)
    
    print('F = ')
    print(numpy.dot(R,P))    
    
    return(R,P)