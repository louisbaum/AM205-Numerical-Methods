# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 13:50:13 2016

@author: Louis
"""

import numpy
import scipy
import matplotlib.pyplot as plt

"""
Plots ||b|| =1 and ||Ab||=1
   parametrizes x,y in terms of theta and r. intersection points have been previously calculated by hand and input into script
"""
def prob1a():
    A = [[3, -1],
         [1, 0]]
    theta = numpy.linspace(0,2*numpy.pi,1000)
    b1x = numpy.cos(theta)
    b1y = numpy.sin(theta)
    r =  1./numpy.sqrt((3*numpy.cos(theta)-numpy.sin(theta))**2 +numpy.cos(theta)**2)#in hindsight I could have used norm(a,2) here
    Ab = [numpy.multiply(r,numpy.cos(theta)),numpy.multiply(r,numpy.sin(theta))]
    fig = plt.gcf()
    plt.plot(b1x,b1y,'red',Ab[0],Ab[1],'black')
    plt.scatter([0,0,2/numpy.sqrt(13),-2/numpy.sqrt(13)],[1,-1,3/numpy.sqrt(13),-3/numpy.sqrt(13)],s=10)
    plt.ylim(-5,5)
    plt.xlim(-5,5)
    plt.title('||b|| shown in red and ||Ab|| shown in black. \n the points b are the intersections ')
    plt.grid()
   
"""
Plots ||c|| =1 and ||Ac||=1
   parametrizes x,y in terms of theta and r. intersection points have been previously calculated by hand and input into script   
"""
def prob1b():
    A = [[3, -1],
         [1, 0]]
    theta = numpy.linspace(0,2*numpy.pi,1000)
    cos = numpy.cos(theta)
    sin = numpy.sin(theta)
    r = 1./ numpy.maximum(numpy.abs(cos),numpy.abs(sin)) #in hindsight I could have used norm(a,inf) here
    x = numpy.multiply(r,cos)
    y = numpy.multiply(r,sin)
                
    Abx = 3*cos-sin
    Aby = cos
        
    rAb =1/ numpy.maximum(numpy.abs(Abx),numpy.abs(Aby))
    Abvalues = [numpy.multiply(rAb,cos),numpy.multiply(rAb,sin)]
    
    fig = plt.gcf()
    plt.plot(x,y,'red',Abvalues[0],Abvalues[1],'black')
    plt.scatter([0,0,2./3.,-2./3.],[1,-1,1,-1],s=10)
    plt.xlim(-5,5)
    plt.ylim(-5,5)
    plt.title('||c||inf shown in red and ||Ac||inf shown in black. \n the points c are the intersections ')
    plt.grid()


"""
Accepts "error" = (acceptible error absolute amount) .001 used to determine intersection points
Performs Newtons method to determine the intersection points of ||d||_4 = 1 = ||Ad||_4
plots ||d||_4 = 1 and  ||Ad||_4=1 emphasizing the intersection points

in brief Newtons method identifies zeroes of a function starting from a guess and using properties of the derivative to iterate towards and answer

F(x+dx) ~ F(x) + F'x *dx

F'(x) is the jacobian (J) for vector systems

but we can flip this around

x+dx ~ x - J^-1 *f(x) or
xi+1 ~ xi - J^-1* F(xi)

we can do this until F(xi) is close to zero (the point we are trying to get too)

"""  
def prob1c(error):
    points =[]
    initialguess = [[2./3.,1.],[0,1.]]
    for i in range(2): #this for loop begine newtons method for two different starting points contained in initialguess.
        guess = initialguess[i] 
        diff = numpy.dot(numpy.abs(function1c(guess)),[1.,1.])
        while (diff > error): #continues iterating until the guess is "close enough" to correct
            fnew = function1c(guess) # calls f(xi) and collects the value
            Jacob= numpy.linalg.inv(Jacobian(guess[0],guess[1])) # calculaes the Jacobian and inverts it
            newguess = numpy.subtract(guess,numpy.dot(Jacob,fnew)) #calculates new guess
            guess = newguess
            diff = numpy.dot(numpy.abs(function1c(guess)),[1,1])
        #once the point is "close enough" stores it and by the symmetry of the problem (the point opposite it)
        points.append(guess)
        points.append(numpy.multiply(guess,[-1,-1]))
        print(guess)
        print(numpy.multiply(guess,[-1,-1]))
        
    #generates the plot    
    theta = numpy.linspace(0,2*numpy.pi,1000)
    r = (numpy.cos(theta)**4+numpy.sin(theta)**4)**(-1./4) #in hindsight I could have used norm(a,4) here
    Ar = ((3*numpy.cos(theta)-numpy.sin(theta))**4+numpy.cos(theta)**4)**(-1./4) #in hindsight I could have used norm(a,4) here
    
    print(numpy.transpose(points))
    plt.plot(numpy.multiply(r,numpy.cos(theta)),numpy.multiply(r,numpy.sin(theta)),'red',numpy.multiply(Ar,numpy.cos(theta)),numpy.multiply(Ar,numpy.sin(theta)),'black')
    plt.scatter(numpy.transpose(points)[0],numpy.transpose(points)[1],s=10)
    plt.ylim(-5,5)
    plt.xlim(-5,5)
    plt.title('||d||4 shown in red and ||Ad||4 shown in black. \n the points d are the intersections ')
    plt.grid()

    
    
    
"""        
Plots the points determined throughout this probelm and draws the lines that intersect them.
"""
def prob1d():
    no1 = [[2/numpy.sqrt(13),-2/numpy.sqrt(13),2./3.,-2./3.,0.63728994,-0.63728994],[3/numpy.sqrt(13),-3/numpy.sqrt(13),1,-1,0.95593491, -0.95593491]]
    ones = [[0,0],[1,-1]]
    plt.scatter(no1[0],no1[1])
    plt.scatter(ones[0],ones[1])
    plt.plot([-2./3,2./3],[-1,1])
    plt.plot([0,0],[1,-1])
    plt.ylim(-1.5,1.5)
    plt.xlim(-1,1)
    plt.title('All points shown. lines x=0 and y=3/2x')
    plt.grid()


"""        
#Accepts a point "d" and calculates the 4 norm for "d" and the 4 norm of "Ad"
#returns a vector of the norm -1 for "d" and "Ad"
"""
def function1c(d):
    Ad =numpy.dot([[3., -1.],[1., 0]],d)    
    nor = numpy.linalg.norm(d,4)
    A4nor = numpy.linalg.norm(Ad,4)
    return([nor-1,A4nor-1])    

"""
Accepts components of a point (x,y) and returns the Jacobian at (x,y).
"""
def Jacobian(x,y):
    xn = 1.*x
    yn = 1.*y
    B=numpy.zeros((2,2))
    B[0][0] = xn**3/(xn**4+yn**4)**(3./4)
    B[0][1] = yn**3/(xn**4+yn**4)**(3./4)
    B[1][0] = (4.*xn**3 + 12*(3*xn-yn)**3)/(4.*(xn**4+(3*xn-yn)**4)**(3./4))
    B[1][1] = -((3.*xn-yn)**3)/((xn**4+(3*xn-yn)**4)**(3./4))
    return((B))
    
    