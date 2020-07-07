# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 13:34:01 2018

@author: kmiri_000
"""
import numpy as np, subs as su, matplotlib.pyplot as plt

#initializing variables
c=2.0/3
n=10
dx=1.0/(n-1)
dt=c*(dx**2)
tmax=1.0
k=1.0
m=int((tmax/dt)+1)

#running the crank scheme 
crank=su.crankn(dx,dt,m,n)

#analytic solution 
ana=su.analyticsol(k,dt,n,dx,tmax)

#defining x range
x=np.arange(0,1.001,dx)

#plotting numerical solution against the analytic one
plt.plot(x,crank[-1],'bx',label="crank")
plt.plot(x,crank[-1],'b')
plt.plot(x,ana[-1],'rx',label="analytic")
plt.plot(x,ana[-1],'r')  
plt.xlabel('x values')
plt.ylabel(r'$\theta$'" values")
plt.legend()
plt.show()

