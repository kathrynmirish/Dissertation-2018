# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 13:28:58 2017

@author: kmiri_000
"""
import subroutines as su
import numpy as np
import matplotlib.pyplot as plt

#varying n and keeping dt fixed
n=4
dt=0.0001
dx=1.0/(n-1)
c=dt/(dx**2)
tmax=0.5
k=1.0

#new arrays
narray=np.array([])
dxarray=np.array([])
msearray2=np.array([])

#varying n values
while n<21:
    
    #using subroutine to calculate root-mse
    mse=su.mse(c,n,dt,dx,tmax,k,0)
    
    #append to array 
    msearray2=np.append(msearray2,mse)
    
    #append n
    narray=np.append(narray,n)
    
    #append dx
    dxarray=np.append(dxarray,dx)
    
    n=n+1
    dx=1.0/(n-1)
    c=dt/(dx**2)
 
#plotting graph of n vs error
plt.plot(narray,msearray2,'bx')
plt.plot(narray,msearray2)
plt.xlabel('N values')
plt.ylabel('E')
plt.show()

#plotting log plot 
plt.plot(np.log(dxarray),np.log(msearray2),'bx')
plt.plot(np.log(dxarray),np.log(msearray2))
plt.xlabel('ln dx values')
plt.ylabel('ln E')
plt.show()

print (np.log(msearray2[-1])-np.log(msearray2[0]))/(np.log(dxarray[-1])-np.log(dxarray[0]))
