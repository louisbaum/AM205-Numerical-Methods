# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 17:39:01 2016

@author: Louis
"""

nr=15
nz=30
nt = 300
tfin = 14.11 #seconds
zfin = 7.62 #meters
T0=273.15
OD = .375*2.54/100 #m
ID = 5./16 *2.54/100 #m
pipe = numpy.zeros([nr+2,nz+2,nt])
rad = numpy.linspace(0,OD/2,nr+2)
copperbound=0
for i in range(nr)[::-1]:
    if(rad[i] < ID/2.):
        copperbound = i+1
        break
times = numpy.linspace(0,tfin,nt)
distance = numpy.linspace(0,zfin,nz)
dz = distance[1]#meters
dr = rad[1] #meters
dt =times[1] # seconds
M = numpy.zeros([(nr+2)*(nz+2)+1,(nr+2)*(nz+2)+1])
vm =.54 # m/s 

