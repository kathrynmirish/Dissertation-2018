# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 13:28:58 2017

@author: kmiri_000
"""
import subroutines as su
import numpy as np
import matplotlib.pyplot as plt

#initialising the variables
c=2.0/3
n=8
dx=1.0/(n-1)
dt=c*dx**2
tmax=1.0
k=1.0

#using the sub routine to get theta matrices
num=su.numericalsol(c,dt,n,tmax)
ana=su.analyticsol(k,dt,n,dx,tmax)

#defining the range for x
x=np.arange(0,tmax+dx/2,dx)

#plotting the graphs for tmax=1
plt.plot(x,ana[-1],'rx',label="analytical")
plt.plot(x,ana[-1],'r')
plt.plot(x,num[-1],'g+',label="numerical")
plt.plot(x,num[-1],'g')
plt.xlabel("x values")
plt.ylabel(r'$\theta$'" values")
plt.legend()
plt.show()

absl=np.abs(ana[-1]-num[-1])
print np.mean(absl)
