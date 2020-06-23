# -*- coding: utf-8 -*-
"""
Created on Fri Dec 09 10:04:31 2016

@author: Louis
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Dec 05 14:54:02 2016

@author: Louis (dicklord69)
"""

import numpy
import scipy.integrate
import matplotlib.pyplot as plt
import matplotlib.animation as an
import time
import random

"""This initialized the functions that return the physical constants as a function of temperature (Kelvin).

thermal conductivity in W/(K-m)
specific heat in J/kg
density kg/m^3
dynamic viscosity in Ns/m^2
kinematic viscosity in m^2/s
"""
execfile("physcon.py")


""" This plot requires running the calculation of imchillcross for nt = 60, 600 and 6000 pipe, beep and boop respectivly.
"""
def comparetimestep():
    plt.figure()
    plt.plot(numpy.linspace(0,60,60),pipe[0,:]-273.15,'red',numpy.linspace(0,60,600),beep[0,:]-273.15,'blue',numpy.linspace(0,60,6000),boop[0,:]-273.15,'green')
    plt.plot(numpy.linspace(0,60,60),pipe[85,:]-273.15,'red',numpy.linspace(0,60,600),beep[85,:]-273.15,'blue',numpy.linspace(0,60,6000),boop[85,:]-273.15,'green')
    plt.plot(numpy.linspace(0,60,60),pipe[160,:]-273.15,'red',numpy.linspace(0,60,600),beep[160,:]-273.15,'blue',numpy.linspace(0,60,6000),boop[160,:]-273.15,'green')
    plt.xlabel("Time (s)")
    plt.ylabel("Temp (C)")
    plt.ylim(0,100)
    plt.legend(["dt = 1s","dt = 0.1 s","dt = 0.01 s"],loc =4)


