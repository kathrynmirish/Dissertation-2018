# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 19:14:10 2018

@author: kmiri_000
"""

def crankn(c,dt,dx,tmax):
    
    import numpy as np

    #defining the constant r
    r=c*dt/dx
    
    #defining n, m
    n=int((4/dx)+1)
    m=int((tmax/dt)+1)
    
    #defining A
    A=np.zeros((n,n))
    
    for i in range(0,n):
        A[i,i]=1
    
    for i in range (0,n-1):
        A[i,i+1]=-r/4
        A[i+1,i]=r/4
    
    #defining A1
    A1=np.matrix.transpose(A)
    
    #inverting A1
    A1inv=np.linalg.inv(A1)
    
    #defining b
    b=np.zeros(n)
    b[0]=-r/2
    b[-1]=r/2

    #defining the u matrix
    u=np.zeros(shape=(m,n))

    #inputting the initial values for u 
    for j in range(0,n):
        u[0,j]=np.sin(2*(-2.0+j*dx))

    #finding the next set of u values by matrix multiplication
    for i in range(0,m-1):
        um=u[i,0:n]
        
        mid=np.dot(A,um)
        
        uplus1=np.dot(A1inv,mid+b)
    
        u[i+1,0:n]=uplus1
    
    return u

