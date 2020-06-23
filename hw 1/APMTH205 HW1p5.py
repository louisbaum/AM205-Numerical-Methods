# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 18:17:50 2016

@author: Louis
"""
import numpy
import scipy
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import PIL

global lsidebyside
global T

#read in images
def a():
    global lsidebyside
    global T
    rightr =scipy.misc.imread(r'right\low3.png')
    rightb =scipy.misc.imread(r'right\low2.png')
    rightg =scipy.misc.imread(r'right\low1.png')
    rreg = scipy.misc.imread(r'right\regular.png')
    
    A= numpy.zeros(((3*len(rreg)*len(rreg[0])),30))
    
    for i in range(len(rreg)):
        for k in range(len(rreg[0])):
            for j in range(3):
                A[3*k+j+i*3*len(rreg[0])][(3*j):(3*(j+1))] = rightr[i][k][0:3] 
                A[3*k+j+i*3*len(rreg[0])][(3*j+9):(3*(j+1)+9)] = rightb[i][k][0:3]
                A[3*k+j+i*3*len(rreg[0])][(3*j+18):(3*(j+1)+18)]= rightg[i][k][0:3]
                A[3*k+j+i*3*len(rreg[0])][27+j] = 1
                
    b = numpy.zeros((3*len(rreg)*len(rreg[0])))
    
    for i in range(len(rreg)):
        for k in range(len(rreg[0])):
            for j in range(3):
                b[j+3*k+i*3*len(rreg[0])] = rreg[i][k][j]
    
    x = numpy.linalg.lstsq(A,b)[0]
    
    fixedpic = numpy.zeros((len(rreg),len(rreg[0]),3))
    for i in range(len(rreg)):
        for k in range(len(rreg[0])):
            for j in range(3):
                fixedpic[i][k][j] = numpy.dot(A[3*k+j+i*3*len(rreg[0])],x)
    fixedpic= numpy.clip(fixedpic,0,255)
    sidebyside = numpy.concatenate((rreg, fixedpic),axis=1)
    
    plt.imshow(sidebyside.astype('uint8'))
    plt.title('original pic on left. reconstructed on right')
    
    S=(numpy.linalg.norm(numpy.subtract(fixedpic,rreg))**2)/(len(rreg)*len(rreg[0]))
    print('S = ')
    print(S)
    
    leftr =scipy.misc.imread(r'left\low3.png')
    leftb =scipy.misc.imread(r'left\low2.png')
    leftg =scipy.misc.imread(r'left\low1.png')
    lreg = scipy.misc.imread(r'left\regular.png')
    
    B = numpy.zeros(((3*len(rreg)*len(rreg[0])),30))
    for i in range(len(rreg)):
        for k in range(len(rreg[0])):
            for j in range(3):
                B[3*k+j+i*3*len(rreg[0])][(3*j):(3*(j+1))] = leftr[i][k][0:3] 
                B[3*k+j+i*3*len(rreg[0])][(3*j+9):(3*(j+1)+9)] = leftb[i][k][0:3]
                B[3*k+j+i*3*len(rreg[0])][(3*j+18):(3*(j+1)+18)]= leftg[i][k][0:3]
                B[3*k+j+i*3*len(rreg[0])][27+j] = 1    
    
    lfixedpic = numpy.zeros((len(lreg),len(lreg[0]),3))
    for i in range(len(lreg)):
        for k in range(len(lreg[0])):
            for j in range(3):
                lfixedpic[i][k][j] = numpy.dot(B[3*k+j+i*3*len(lreg[0])],x)
    lfixedpic= numpy.clip(lfixedpic,0,255)
    lsidebyside = numpy.concatenate((lreg, lfixedpic),axis=1)
    T =(numpy.linalg.norm(numpy.subtract(lfixedpic,lreg))**2)/(len(lreg)*len(lreg[0]))
    
def b(): 
    global lsidebyside
    global T
    plt.imshow(lsidebyside.astype('uint8'))
    plt.title('original pic on left. reconstructed on right')
    print('T = ')
    print(T)
    