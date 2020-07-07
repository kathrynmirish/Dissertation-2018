# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 10:16:12 2018

@author: kmiri_000
"""
import subroutines as su
import numpy as np
import matplotlib.pyplot as plt

n=21
dx=1.0/(n-1)
tmax=0.5
k=1.0

#defining arrays for graph
tarray=np.array([1.0/850,1.0/840,1.0/830,1.0/820])
msearray=np.array([])

#varying c values
for i in range(0,4):
    
    c=tarray[i]/(dx**2)
    
    #using subroutine to calculate root mean squared error
    mse=su.mse(c,n,tarray[i],dx,tmax,k,0)
    
    #append to the array 
    msearray=np.append(msearray,mse)
    


print (np.log(msearray[-1])-np.log(msearray[1]))/(np.log(tarray[-1])-np.log(tarray[1]))

#plotting graph of c vs error    
plt.plot(tarray,msearray)
plt.plot(tarray,msearray,'bx')
plt.xlabel("dt values")
plt.ylabel("E")
plt.show()


#rounding error example
n=21
dx=1.0/(n-1)
tmax=0.5
k=1.0

#defining arrays for graph
tarray=np.array([0.00001,0.00005,0.0001])
msearray=np.array([])

#varying c values
for i in range(0,3):
    
    c=tarray[i]/(dx**2)
    
    #using subroutine to calculate root mean squared error
    mse=su.mse(c,n,tarray[i],dx,tmax,k,0)
    
    #append to the array 
    msearray=np.append(msearray,mse)
    


print (msearray[-1]-msearray[1])/(tarray[-1]-tarray[1])

#plotting graph of c vs error    
plt.plot(tarray,msearray)
plt.plot(tarray,msearray,'bx')
plt.xlabel("dt values")
plt.ylabel("E")
plt.show()
