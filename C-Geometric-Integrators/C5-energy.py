# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 11:32:16 2018

@author: kmiri_000
"""

import numpy as np, subs as su, matplotlib.pyplot as plt

#FORWARD EULER
#intializing variables
e=0.9
p=2*np.pi
N=100
h=p/300

#defining initial conditions for t=0
r=[1.0+e,0.0]
v=[0.0,np.sqrt((1.0-e)/(1.0+e))]

#defining empy arrays for numerical solution
rarray=np.array(r)
varray=np.array(v)

#number of steps to make a n orbits
n=int((N*p)/h)+1

#looping to determine n new points
for i in range(0,n-1):
    
    #using the time step to move on one step
    newr=su.euler(h,rarray[-2:],varray[-2:])[0:2]
    newv=su.euler(h,rarray[-2:],varray[-2:])[2:4]
    
    #addinf the new points to the array
    rarray=np.append(rarray,newr)
    varray=np.append(varray,newv)
    
#defining array of time
tarray=np.arange(0,(n-1)*h,h)

#defining E0
E0=0.5*((1.0-e)/(1.0+e))-(1/np.linalg.norm(r))

#defining fractional error in energy for forward euler
Eerror=np.array([])

for i in range(0,n-1):
    Ei=0.5*((varray[2*i])**2+(varray[2*i+1])**2)-(1/np.linalg.norm([rarray[2*i],rarray[2*i+1]]))
    
    frac=(np.absolute(Ei-E0))/np.absolute(E0)
    
    Eerror=np.append(Eerror,frac)

#MODIFIED EULER
#intializing variables
e=0.9
p=2*np.pi
N=100
h=p/150

#defining initial conditions for t=0
r=[1.0+e,0.0]
v=[0.0,np.sqrt((1.0-e)/(1.0+e))]

#defining empy arrays for numerical solution
rarray2=np.array(r)
varray2=np.array(v)

#number of steps to make a n orbits
n=int((N*p)/h)+1

#looping to determine n new points
for i in range(0,n-1):
    
    #using the time step to move on one step
    newr=su.modeuler(h,rarray2[-2:],varray2[-2:])[0:2]
    newv=su.modeuler(h,rarray2[-2:],varray2[-2:])[2:4]
    
    #addinf the new points to the array
    rarray2=np.append(rarray2,newr)
    varray2=np.append(varray2,newv)

#defining fractional error in energy for mod euler
Eerror2=np.array([])

for i in range(0,n-1):
    Ei=0.5*((varray2[2*i])**2+(varray2[2*i+1])**2)-(1/np.linalg.norm([rarray2[2*i],rarray2[2*i+1]]))
    
    frac=(np.absolute(Ei-E0))/np.absolute(E0)
    
    Eerror2=np.append(Eerror2,frac)
    
#LEAPFROG
#intializing variables
e=0.9
p=2*np.pi
N=100
h=p/150

#defining initial conditions for t=0
r=[1.0+e,0.0]
v=[0.0,np.sqrt((1.0-e)/(1.0+e))]

#defining empy arrays for numerical solution
rarray3=np.array(r)
varray3=np.array(v)

#number of steps to make a n orbits
n=int((N*p)/h)+1

#looping to determine n new points
for i in range(0,n-1):
    
    #using the time step to move on one step
    newr=su.leapfrog(h,rarray3[-2:],varray3[-2:])[0:2]
    newv=su.leapfrog(h,rarray3[-2:],varray3[-2:])[2:4]
    
    #addinf the new points to the array
    rarray3=np.append(rarray3,newr)
    varray3=np.append(varray3,newv)

#defining array of time
tarray2=np.arange(0,(n-1)*h,h)

#defining fractional error in energy for mod euler
Eerror3=np.array([])

for i in range(0,n-1):
    Ei=0.5*((varray3[2*i])**2+(varray3[2*i+1])**2)-(1/np.linalg.norm([rarray3[2*i],rarray3[2*i+1]]))
    
    frac=(np.absolute(Ei-E0))/np.absolute(E0)
    
    Eerror3=np.append(Eerror3,frac)
    
#RK4
#intializing variables
e=0.9
p=2*np.pi
N=100
h=p/75

#defining initial conditions for t=0
r=[1.0+e,0.0]
v=[0.0,np.sqrt((1.0-e)/(1.0+e))]

#defining empy arrays for numerical solution
rarray4=np.array(r)
varray4=np.array(v)

#number of steps to make a n orbits
n=int((N*p)/h)+1

#looping to determine n new points
for i in range(0,n-1):
    
    #using the time step to move on one step
    newr=su.RK4(h,rarray4[-2:],varray4[-2:])[0:2]
    newv=su.RK4(h,rarray4[-2:],varray4[-2:])[2:4]
    
    #addinf the new points to the array
    rarray4=np.append(rarray4,newr)
    varray4=np.append(varray4,newv)

#defining array of time
tarray3=np.arange(0,(n-1)*h,h)

#defining fractional error in energy for mod euler
Eerror4=np.array([])

for i in range(0,n-1):
    Ei=0.5*((varray4[2*i])**2+(varray4[2*i+1])**2)-(1/np.linalg.norm([rarray4[2*i],rarray4[2*i+1]]))
    
    frac=(np.absolute(Ei-E0))/np.absolute(E0)
    
    Eerror4=np.append(Eerror4,frac)
    
plt.plot(tarray[0:3000:200],Eerror[0:3000:200],'bx',label="Forward Euler")
plt.plot(tarray[0:3000:200],Eerror[0:3000:200],'b')
plt.plot(tarray2[0:1500:100],Eerror2[0:1500:100],'rx',label="Modified Euler")
plt.plot(tarray2[0:1500:100],Eerror2[0:1500:100],'r')
plt.plot(tarray2[0:1500:100],Eerror3[0:1500:100],'gx',label="Leapfrog")
plt.plot(tarray2[0:1500:100],Eerror3[0:1500:100],'g')
plt.plot(tarray3[0:750:50],Eerror4[0:750:50],'cx',label="RK4")
plt.plot(tarray3[0:750:50],Eerror4[0:750:50],'c')
plt.xlabel("t")
plt.ylabel("fractional error in energy")
plt.legend()
plt.show()
