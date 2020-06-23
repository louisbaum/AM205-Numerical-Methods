# -*- coding: utf-8 -*-
"""
Created on Tue Dec 06 14:53:52 2016

@author: Louis
"""

import numpy
import scipy
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import matplotlib.animation as an
import time
import random

"""Data
"""

#Temp(K), density(kg/m^3),specific heat(J/kg-K), dynamic viscosity (Ns/m^2)
watertable =[[0+273.15,5+273.15,10+273.15,20+273.15,30+273.15,40+273.15,50+273.15,60+273.15,70+273.15,80+273.15,90+273.15,100+273.15], [999.8,999.9,999.7,998.2,995.7,992.2,988.1,983.2,977.8,971.8,965.3,958.4],[4217,4202,4194,4182,4178,4179,4182,4185,4191,4198,4208,4219],[1.792e-3,1.519e-3,1.307e-3,1.002e-3,.798e-3,.653e-3,.547e-3,.467e-3,.404e-3,.355e-3,.315e-3,.282e-3]]
#Temp(K), thermal conductivity (W/mK)
watercptable=[[275.00,280.00,285.,290.,295.,300.,305.,310.,315.,320.,325.,330.,335.,340.,345.,350.,355.,360.,365.,370.],[.5606,.5715,.5818,.5917,.6009,.6096,.6176,.6252,.6322,.6387,.6445,.6499,.6546,.6588,.6624,.6655,.6680,.6700,.6714,.6723]]

#temp(K), cp(J/kg-K)
coppercptable = [[260,280,300,400],[.376e3,.381e3,.386e3,.396e3]]

"""Creating the interpolation functions
"""
waterk =scipy.interpolate.CubicSpline(watercptable[0],watercptable[1])
watercp= scipy.interpolate.UnivariateSpline(watertable[0],watertable[2],s=1)
coppercp = scipy.interpolate.CubicSpline(coppercptable[0],coppercptable[1])
waterp = scipy.interpolate.CubicSpline(watertable[0],watertable[1])
watermu= scipy.interpolate.CubicSpline(watertable[0],watertable[3])

"""return k, the thermal conductivity in W/(K-m)
""" 
def k(T):
    water =float(waterk(T))
    copper=400. #W/(K-m)
    stainless = 16.5
    wort = water*1.060
    return(water,wort,copper,stainless)

"""return cp, the specific heat in J/kg
"""
def cp(T):
    water = float(watercp(T))
    copper =float(coppercp(T))
    wort = float(water)
    stainless =520
    return(water,wort,copper,stainless)
    
"""return p, the density kg/m^3
"""
def p(T):
    water = float(waterp(T))
    copper =8940
    wort = water*1.060
    stainless =7750
    return(water,wort,copper,stainless)
    
"""return the dynamic viscosity in Ns/m^2
"""
def mu(T):
    water= float(watermu(T))
    wort=1.7*water
    return(water,wort)
    
"""return kinematic viscosity in m^2/s
"""
def viscosity(T):
    water,wort = mu(T)
    pwater,pwort,dum1,dum2 = p(T)
    return(water/pwater,wort/pwort)
    
    
"""Plots the interpolated functions to make sure it looks right
"""
    
def checkvalues():
    testk = numpy.zeros([300,4])
    testcp=numpy.zeros([300,4])
    testp=numpy.zeros([300,4])
    testmu=numpy.zeros([300,2])
    testvis = numpy.zeros([300,2])
    Tvals  = numpy.linspace(273.15,373,300)
    for i in range(300):
        testk[i,:]=k(Tvals[i])
        testcp[i,:] = cp(Tvals[i])
        testp[i,:] = p(Tvals[i])
        testmu[i,:] = mu(Tvals[i])
        testvis[i,:] = viscosity(Tvals[i])
    for i in range(4):
        plt.figure()
        plt.plot(Tvals, testcp[:,i])
        plt.ylabel('testcp')
    for i in range(4):
        plt.figure()
        plt.plot(Tvals, testk[:,i])
        plt.ylabel('testk')
    for i in range(4):
        plt.figure()
        plt.plot(Tvals, testp[:,i])
        plt.ylabel('testp')
    for i in range(2):
        plt.figure()
        plt.plot(Tvals, testmu[:,i])
        plt.ylabel('testmu')
    for i in range(2):
        plt.figure()
        plt.plot(Tvals, testvis[:,i])
        plt.ylabel('viscosity')
 