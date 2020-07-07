
"""
Created on Fri Jan 26 19:19:44 2018

@author: kmiri_000
"""

import numpy as np, subroutines as su, matplotlib.pyplot as plt

#intializing variables
n=21
dx=1.0/(n-1)
tmax=0.5
k=1.0

#defining arrays for graph
tarray=np.array([1.0/30,1.0/28,1.0/24,1.0/20,1.0/16,1.0/10])
msearray1=np.array([])

#varying c values
for i in range(0,6):
    
    #define c
    c=((n-1)**2)*tarray[i]
    
    #using subroutine to calculate root mean squared error
    mse1=su.mse(c,n,tarray[i],dx,tmax,k,1)
    
    #append to the array 
    msearray1=np.append(msearray1,mse1)


#log plot 
plt.plot(np.log(tarray),np.log(msearray1),'bx',label='crank')
plt.plot(np.log(tarray),np.log(msearray1),'b')
plt.xlabel('dt values')
plt.ylabel('E')
plt.legend()
plt.show()

print (np.log(msearray1[-1])-np.log(msearray1[0]))/(np.log(tarray[-1]-np.log(tarray[0])))

#intializing variables
dt=0.00005
tmax=0.5
k=1.0

#defining arrays for graph
dxarray=np.array([1.0/5,1.0/6,1.0/7,1.0/8,1.0/9,1.0/10])
narray=np.array([6,7,8,9,10,11])
msearray1=np.array([])

#varying c values
for i in range(0,6):
    
    #using subroutine to calculate root mean squared error
    mse1=su.mse(c,narray[i],dt,dxarray[i],tmax,k,1)
    
    #append to the array 
    msearray1=np.append(msearray1,mse1)
    

print (np.log(msearray1[-1])-np.log(msearray1[0]))/(np.log(dxarray[-1])-np.log(dxarray[0]))
