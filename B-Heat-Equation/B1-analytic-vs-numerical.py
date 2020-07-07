# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 13:28:58 2017

@author: kmiri_000
"""
import subroutines as su
import numpy as np
import matplotlib.pyplot as plt

#initialising the variables
c=1.0/3
n=5
dx=1.0/(n-1)
dt=c*dx**2
tmax=1.0
k=1.0

#using the sub routine to get theta matrices
num=su.numericalsol(c,dt,n,tmax)
ana=su.analyticsol(k,dt,n,dx,tmax)

print num[36], ana[36], num[36]-ana[36]

#defining the range for x
x=np.arange(0,tmax+dx/2,dx)

#plotting the graphs 12*dt=0.25
plt.plot(x,ana[12],'rx',label="analytical 0.25")
plt.plot(x,num[12],'g+',label="numerical 0.25")
plt.xlabel("x values")
plt.ylabel(r'$\theta$'" values")
plt.plot(x,ana[36],'bx',label="analytical 0.75")
plt.plot(x,num[36],'c+',label="numerical 0.75")
plt.legend()
plt.show()