# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 19:28:39 2018

@author: kmiri_000
"""

import wavesubs as su, numpy as np, matplotlib.pyplot as plt

c=0.5
tmax=0.3

crank=su.crankn(c,0.001,0.25,tmax)

print crank[-1]

x=np.arange(-2,2.125,0.25)

analytic=np.array([])

for i in range(0,len(x)):
    ana=np.sin(2*x[i]+c*tmax)
    analytic=np.append(analytic,ana)
    
print analytic

plt.plot(x,analytic,'bx',label="explicit")
plt.plot(x,analytic,'b')
plt.plot(x,crank[-1],'rx',label="crank")
plt.plot(x,crank[-1],'r')
plt.xlabel("x")
plt.ylabel("u")
plt.legend()
plt.show()
