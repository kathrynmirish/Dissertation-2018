#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 12:19:48 2017

@author: amtsdg
"""

import numpy as np, subroutines as su, scipy, matplotlib.pyplot as plt


# parameters
dt=0.01      # time step
a=np.pi/2    # initial angle
nonlin=1      # 0 or 1, for linear or nonlinear

dtr=np.array([])
errorr=np.array([])

while dt<0.5:
    
    # initialising variables
    t=0                    # start time at 0
    theta=np.array([a,0])  # theta[0]=theta, theta[1]=\dot theta
    n=0    
    i=0                    # number of steps taken

    # record t and theta
    tr=np.array([])
    thetar=np.array([])

    #repeats until the sign changes
    while i<2:
        
        #rk3 for other graph
        theta=su.step_rk2a(theta,dt,nonlin)
            
        t=t+dt
        
        n=n+1
        
        tr=np.append(tr,t)
        thetar=np.append(thetar,theta[0])
        
        if theta[0]<0:
            i=i+1
    
    #linear interpolation of final 2 points to get estimate of period
    tapprox = su.linin(tr,thetar)
    
    #exact period value for a=pi/2
    #exactt=2*np.pi
    
    #exact period value non-linear
    print tapprox
    
    exactt = (scipy.special.gamma(0.25))**2/(np.sqrt(np.pi))
    
    #error 
    error = abs(tapprox - exactt)

    dtr=np.append(dtr,dt)
    errorr=np.append(errorr,error)
    
    dt=dt*2

#print dtr, errorr

logdt=np.log2(dtr)
logerror=np.log2(errorr)

#convergence plot 
plt.plot(logdt,logerror,'bx')
plt.plot(logdt,logerror)
plt.title("")
plt.xlabel("value of " r'$\log_2\Delta t$')
plt.ylabel("value of " r'$\log_2 E$')
plt.show()

print (logerror[-1]-logerror[1])/(logdt[-1]-logdt[1])
    