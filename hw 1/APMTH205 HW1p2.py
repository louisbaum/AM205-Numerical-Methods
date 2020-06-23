# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 10:15:13 2016

@author: Louis
"""
import numpy
import scipy
import matplotlib.pyplot as plt

points =[]
n=4
numsamples = 1000

def f(x):
    return(numpy.exp(4*x)+numpy.exp(-2*x))

def cheby(n):
    for j in range(1,n+1):
        xj=numpy.cos((2*j-1)*numpy.pi/(2*n))
        points.append([xj,f(xj)])
    
def lagrange(x):
    cheby(n)
    poly = [[0,0,0,0]]*n
    prodvalue =[0]*n
    for k in range(n):
        poly[k] = [(x-points[j][0])/(points[k][0]-points[j][0]) for j in range(n) if j !=k]
    for i in range(n):
        prodvalue[i] = points[i][1]*numpy.prod(poly[i])
    return([numpy.sum(prodvalue),prodvalue[1]])

xvalues = numpy.linspace(-1,1,numsamples)
fyvalues = f(xvalues)
lagrangevalues = []
for i in range(numsamples):
    lagrangevalues.append(lagrange(xvalues[i]))
    

plt.plot(xvalues,fyvalues,'blue',xvalues,numpy.transpose(lagrangevalues)[0],'red')
plt.title("f(x) in blue, lagrange polynomail interpolation in red")

#plt.plot(xvalues,abs(fyvalues-lagrangevalues))

print("The maximum ||f-p3||(inf) is :")
print(numpy.amax(abs(fyvalues-numpy.transpose(lagrangevalues)[0])))
