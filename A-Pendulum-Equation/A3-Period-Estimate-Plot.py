#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 12:19:48 2017

@author: amtsdg
"""

import numpy as np, subroutines as su, matplotlib.pyplot as plt

# parameters
dt=0.001      # time step
a=(np.pi)/48    # initial angle
nonlin=1      # 0 or 1, for linear or nonlinear
tarray=np.array([])
aarray=np.array([])

while a<(5*np.pi)/6:        
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
            
        theta=su.step_rk4(theta,dt,nonlin)
                
        t=t+dt
            
        n=n+1
            
        tr=np.append(tr,t)
        thetar=np.append(thetar,theta[0])
            
        if theta[0]<0:
            i=i+1
        
    #linear interpolation of final 2 points to get estimate of period
    tapprox = su.cubein(tr,thetar)
    tarray=np.append(tarray,tapprox)
    aarray=np.append(aarray,a/np.pi)
    
    a=a+(np.pi/24)
    
piline = np.array([2*np.pi for j in xrange(len(tarray))])
    
plt.plot(aarray,tarray)
plt.plot(aarray,tarray,'bx',label='T(a)')
plt.plot(aarray,piline,'r--',label='Limit of ' r'$2\pi$')
plt.xlabel("a"r'$/\pi$')
plt.ylabel("estimate of T(a)")
plt.legend()
plt.show

