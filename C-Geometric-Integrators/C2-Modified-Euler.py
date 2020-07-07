# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 15:30:15 2018

@author: kmiri_000
"""

import numpy as np, subs as su, matplotlib.pyplot as plt

#intializing variables
e=0.5
h=(2.0*np.pi)/500
N=1

#defining initial conditions for t=0
r=[1.0+e,0.0]
v=[0.0,np.sqrt((1.0-e)/(1.0+e))]

#defining empy arrays for numerical solution
rarray=np.array(r)
varray=np.array(v)

#number of steps to make a full orbit
n=int((N*2*np.pi)/h)+1

#looping to determine n new points
for i in range(0,n-1):
    
    #using the time step to move on one step
    newr=su.modeuler(h,rarray[-2:],varray[-2:])[0:2]
    newv=su.modeuler(h,rarray[-2:],varray[-2:])[2:4]
    
    #addinf the new points to the array
    rarray=np.append(rarray,newr)
    varray=np.append(varray,newv)

#defining empty array for the analytical solution
rana=np.array(r)

for i in range(1,n):
    
    #defining E
    E=i*h
    
    #returning new values for r
    newr=su.analytical(e,E)
    
    #appending to the array
    rana=np.append(rana,newr)
    
#slicing into x and y values only using a fraction of the points
xnum=rarray[0:len(rarray):8]
ynum=rarray[1:len(rarray):8]
xana=rana[0:len(rana):8]
yana=rana[1:len(rana):8]

#plotting the graph
plt.scatter(xnum,ynum,s=5,label="numerical solution")
plt.scatter(xana,yana,s=5,label="analytic solution")
plt.scatter(0.0,0.0,s=50,label="position of star")
plt.xlabel("x axis")
plt.ylabel("y axis")
plt.legend()
plt.show()

