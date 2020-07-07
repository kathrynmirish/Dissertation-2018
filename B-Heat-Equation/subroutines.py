# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 16:35:41 2017

@author: kmiri_000
"""

def analyticsol(k,dt,n,dx,tmax):
#returns the analytic solution of the equation

    import numpy as np
    
    #define m as number of points in [0,tmax]
    m=int(tmax/dt+1)

    #defining a mxn matrix for the theta values 
    theta=np.zeros(shape=(m,n))
    
    #filling the matrix with theta values 
    for i in range(0,m):
        for j in range (0,n):
            #analytic solution
            theta[i][j]=np.sin(np.pi*j*dx)*np.exp(-k*(np.pi**2)*i*dt)
    
    return theta 


def numericalsol(c,dt,n,tmax):
#finding the numerical solution

    import numpy as np
    
    #define m as number of points in [0,tmax]
    m=int(tmax/dt+1)
    
    #define dx
    dx=1.0/(n-1)
    
    #defining the mxn matrix for the theta values
    theta=np.zeros(shape=(m,n))
    
    #intial conditions for theta=f(x)
    for j in range(1,n-1):
        theta[0][j]=np.sin(j*np.pi*dx)
    
    #using the numerical scheme 
    for i in range(1,m):
            for j in range(1,n-1):
                theta[i][j]=(theta[i-1][j]+c*(theta[i-1][j+1]-2*theta[i-1][j]+theta[i-1][j-1]))
    
    return theta

def crank(dt,dx,m):
    
    import numpy as np

    #defining the constant r
    r=dt/(dx**2)
    r=r/2
    
    #defining D and I
    D=np.matrix('2,-1,0;-1,2,-1;0,-1,2')
    I=np.eye(3)

    #defining the matrix A1
    A1=I+r*D
    A1=np.linalg.inv(A1)

    #defining the matrix A2
    A2=I-r*D

    #defining A
    A=np.dot(A1,A2)

    #defining the theta matrix
    theta=np.zeros(shape=(m,5))

    #inputting the initial values for theta 
    for j in range(1,4):
        theta[0,j]=np.sin(j*np.pi*dx)

    #finding the next set of theta values by matrix multiplication
    for i in range(0,m-1):
        thetam=theta[i,1:4]
        
        thetamplus1=np.dot(A,thetam)
    
        theta[i+1,1:4]=thetamplus1
    
    return theta


def crank10(dt,dx,m):
    
    import numpy as np

    #defining the constant r
    r=(dt/(dx**2))/2

    #defining D and I
    D=np.matrix('2,-1,0,0,0,0,0,0;-1,2,-1,0,0,0,0,0;\
                0,-1,2,-1,0,0,0,0;0,0,-1,2,-1,0,0,0;\
                0,0,0,0,-1,2,-1,0;0,0,0,0,0,-1,2,-1;\
                0,0,0,0,0,0,-1,2;0,0,0,0,0,0,0,-1')
    I=np.eye(8)

    #defining the matrix A1
    A1=np.linalg.inv(I+r*D)

    #defining the matrix A2
    A2=I-r*D

    #defining A
    A=np.dot(A1,A2)

    #defining the theta matrix
    theta=np.zeros(shape=(m,10))

    #inputting the initial values for theta 
    for j in range(1,9):
        theta[0][j]=np.sin(j*np.pi*dx)

    #finding the next set of theta values by matrix multiplication
    for i in range(0,m-1):
        thetam=theta[i][1:9]
        
        thetamplus1=np.dot(A,thetam)
    
        theta[i+1][1:9]=thetamplus1
    
    return theta

def crankn(dt,dx,m,n):
    
    import numpy as np

    #defining the constant r
    r=dt/(dx**2)
    r=r/2
    
    #defining D and I
    D=np.zeros((n-2,n-2))
    
    for i in range(0,n-2):
        D[i,i]=2
    
    for i in range (0,n-3):
        D[i,i+1]=-1
        D[i+1,i]=-1
        
    I=np.eye(n-2)

    #defining the matrix A1
    A1=I+r*D
    A1=np.linalg.inv(A1)

    #defining the matrix A2
    A2=I-r*D

    #defining A
    A=np.dot(A1,A2)

    #defining the theta matrix
    theta=np.zeros(shape=(m,n))

    #inputting the initial values for theta 
    for j in range(1,n-1):
        theta[0,j]=np.sin(j*np.pi*dx)

    #finding the next set of theta values by matrix multiplication
    for i in range(0,m-1):
        thetam=theta[i,1:n-1]
        
        thetamplus1=np.dot(A,thetam)
    
        theta[i+1,1:n-1]=thetamplus1
    
    return theta




def mse(c,n,dt,dx,tmax,k,approx):
#calculates root mean squared error
#approx is 0 for numericalsol, 1 for crank nicholson
    
    if approx==0:
        #numerical solutions matrix
        num=numericalsol(c,dt,n,tmax)
        
    elif approx==1:
        #define m
        m=int((tmax/dt)+1)
        
        #crank matrix
        num=crankn(dt,dx,m,n)
    
    #analytical solutions matrix 
    ana=analyticsol(k,dt,n,dx,tmax)
    
    #error matrix 
    error=abs(num-ana)
    
    #removing the initial values
    error=error[1:,1:n-1]
    
    #mean of error values 
    mse=error.mean()
    
    return mse