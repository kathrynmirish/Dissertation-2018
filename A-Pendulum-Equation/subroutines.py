#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 12:42:27 2017

@author: amtsdg
"""

def theta_dot(theta,nonlin):
    
    import numpy as np
    
    dtheta=0*theta # initialise as array of correct size
    
    dtheta[0]=theta[1]
    if nonlin == 0:
        dtheta[1]=-theta[0]
    else:
        dtheta[1]=-np.sin(theta[0])
        
    return dtheta

def step_euler(theta,dt,nonlin):
        
    # Euler time step for theta[0,1] by time dt
    
    # calculate time derivative
    dtheta=theta_dot(theta,nonlin)
    
    # update
    theta_out=theta+dt*dtheta

    return theta_out

def step_rk2a(theta,dt,nonlin):
        
    # RK2 time step for theta[0,1] by time dt
    
    # calculate time derivative
    dtheta1=theta_dot(theta,nonlin)
    
    # first trial step
    theta1=theta+dt*dtheta1
    
    # calculate time derivative 2
    dtheta2=theta_dot(theta1,nonlin)
    
    # update
    theta_out=theta+0.5*dt*(dtheta1+dtheta2)
    
    return theta_out

def step_rk2b(theta,dt,nonlin):
        
    # RK2 time step for theta[0,1] by time dt
    
    # calculate time derivative
    k1=theta_dot(theta,nonlin)
    
    # calculate time derivative 2
    k2=theta_dot(theta+(2/3)*dt*k1,nonlin)
    
    # update
    theta_out=theta+dt*(0.25*k1+0.75*k2)
    
    return theta_out

#I wrote this part 
def step_rk3(theta,dt,nonlin):
    
    #RK3 time step for theta[0,1] by time dt
    
    k1=theta_dot(theta,nonlin)
    
    k2=theta_dot(theta+k1*dt*0.5,nonlin)

    k3=theta_dot(theta-k1*dt+2*k2*dt,nonlin)
    
    theta_out=theta+dt*(k1+4*k2+k3)/6
    
    return theta_out

#I wrote this part
def step_rk4(theta,dt,nonlin):
    
    #RK4 time step for theta[0,1] by time dt
    
    k1=theta_dot(theta,nonlin)
    
    k2=theta_dot(theta+k1*dt*0.5,nonlin)
    
    k3=theta_dot(theta+k2*dt*0.5,nonlin)
    
    k4=theta_dot(theta+k3*dt,nonlin)
    
    #update
    theta_out=theta+dt*(k1+2*k2+2*k3+k4)/6
    
    return theta_out

def test_rk(dt,rk):

    import numpy as np
    #returns error of different RKs with different time steps
    #starting angle is 1, linear equation only, tmax=5
    #dt time step
    #rk 1, 2, 3, or 4

    # initialising variables
    t=0                    # start time at 0
    theta=np.array([1,0])  # theta[0]=theta, theta[1]=\dot theta
    n=0                    # number of steps taken

    if rk==1:
        timestep=step_euler
    elif rk==2:
        timestep=step_rk2a
    elif rk==3:
        timestep=step_rk3
    elif rk==4:
        timestep=step_rk4
    
    
    nsteps=round(5/dt) # number of steps
    
    while n < nsteps:
        
        theta=timestep(theta,dt,0)
            
        t=t+dt
        
        n=n+1
    
    theta_exact=np.cos(t)
    error=np.abs(theta_exact-theta[0])
    return(error)


#I also wrote this part
def linin(tr,thetar):
    
    #use linear interpolation to work out quarter step
    theta=thetar[-3:-1]
    t=tr[-3:-1]
    
    tquart = t[0] - theta[0]*(t[1]-t[0])/(theta[1]-theta[0])
    
    #return full step approximation
    return tquart*4

def quadin(tr,thetar):
    import numpy as np
    
    #neater version of original code quadin2 
    #lagrange interpolation with 3 points
    
    #slice the arrays to get final 3 points
    t=tr[-4:-1]
    theta=thetar[-4:-1]
    
    #coefficients of the polynomial
    coeff2=0
    coeff1=0
    coeff0=0
    
    for i in range(0,3):
        #iterate to get coefficient for x^2
        coeff2=coeff2+(theta[i]/((t[i]-t[np.mod(i+1,3)])\
        *(t[i]-t[np.mod(i+2,3)])))
        
        #iterate to get coefficient for x
        coeff1=coeff1-theta[i]*(t[np.mod(i+1,3)]+t[np.mod(i+2,3)])\
        /((t[i]-t[np.mod(i+1,3)])*(t[i]-t[np.mod(i+2,3)]))
        
        #iterate to get coefficient for 1
        coeff0=coeff0+(theta[i]*t[np.mod(i+1,3)]*t[np.mod(i+2,3)])\
        /((t[i]-t[np.mod(i+1,3)])*(t[i]-t[np.mod(i+2,3)]))
        
    #find roots of polynomial
    root=np.polynomial.polynomial.polyroots([coeff0,coeff1,coeff2])
    
    #return root in correct range 
    if t[1]<root[0]<t[2]:
        return root[0]*4
    elif t[1]<root[1]<t[2]:
        return root[1]*4

def quadin2(tr,thetar):
    import numpy as np
    
    #lagrange interpolation with 3 points
    
    #picking out the 3 points
    theta1=thetar[-4]
    theta2=thetar[-3]
    theta3=thetar[-2]
    t1=tr[-4]
    t2=tr[-3]
    t3=tr[-2]
    
    #defining the coefficients of the quadratic
    
    #coefficient of t^2
    coeff2=theta1/((t1-t2)*(t1-t3))+\
    theta2/((t2-t1)*(t2-t3))+\
    theta3/((t3-t1)*(t3-t2))
    
    #coefficient of t
    coeff1=(theta1*(-1)*(t2+t3))/((t1-t2)*(t1-t3))+\
    (theta2*(-1)*(t1+t3))/((t2-t1)*(t2-t3))+\
    (theta3*(-1)*(t1+t2))/((t3-t1)*(t3-t2))
    
    #coefficient of 1
    coeff0=(theta1*t2*t3)/((t1-t2)*(t1-t3))+\
    (theta2*t1*t3)/((t2-t1)*(t2-t3))+\
    (theta3*t1*t2)/((t3-t1)*(t3-t2))
    
    
    #find roots of polynomial
    root=np.polynomial.polynomial.polyroots([coeff0,coeff1,coeff2])

    #return root in correct range 
    if t2<root[0]<t3:
        return root[0]*4
    elif t2<root[1]<t3:
        return root[1]*4

    
def cubein(tr,thetar):
    
    import numpy as np
    
    #slice the arrays
    t=tr[-4:]
    theta=thetar[-4:]

    #define 4 variables 
    coeff3=0
    coeff2=0
    coeff1=0
    coeff0=0
    
    #iterate to define the variables
    for i in range(0,4):
        
        coeff3=coeff3+theta[i]/((t[i]-t[np.mod(i+1,4)])\
        *(t[i]-t[np.mod(i+2,4)])*(t[i]-t[np.mod(i+3,4)]))
    
        coeff2=coeff2-theta[i]*(t[np.mod(i+1,4)]+t[np.mod(i+2,4)]\
        +t[np.mod(i+3,4)])/((t[i]-t[np.mod(i+1,4)])\
        *(t[i]-t[np.mod(i+2,4)])*(t[i]-t[np.mod(i+3,4)]))
        
        coeff1=coeff1+theta[i]*(t[np.mod(i+1,4)]*t[np.mod(i+2,4)]\
        +t[np.mod(i+1,4)]*t[np.mod(i+3,4)]+t[np.mod(i+2,4)]\
        *t[np.mod(i+3,4)])/((t[i]-t[np.mod(i+1,4)])\
        *(t[i]-t[np.mod(i+2,4)])*(t[i]-t[np.mod(i+3,4)]))
    
        coeff0=coeff0-theta[i]*((t[np.mod(i+1,4)]*t[np.mod(i+2,4)]\
        *t[np.mod(i+3,4)]))/((t[i]-t[np.mod(i+1,4)])*(t[i]-t[np.mod(i+2,4)])\
        *(t[i]-t[np.mod(i+3,4)]))
        
    #find roots of polynomial
    root=np.polynomial.polynomial.polyroots([coeff0,coeff1,coeff2,coeff3])

    if t[1]<root[0]<t[1]:
       return root[0]*4
    elif t[1]<root[1]<t[2]:
        return root[1]*4
    elif t[1]<root[2]<t[2]:
        return root[2]*4
  
def polyin(tr,thetar,d):
#same as linin/ quadin/ cubein but using inbuilt functions
   
    import numpy as np
    
    #d is degree 1, 2, or 3.
    #slicing the arrays and using polyfit to fit a polynomial
    if d==5:
        theta=thetar[-6:]
        t=tr[-6:]
        pol=np.polyfit(t,theta,5)
    
    if d==4:
        theta=thetar[-5:]
        t=tr[-5:]
        pol=np.polyfit(t,theta,4)
    
    if d==3:
        theta=thetar[-4:]
        t=tr[-4:]
        pol=np.polyfit(t,theta,3)
        
    elif d==2:
        theta=thetar[-4:-1]
        t=tr[-4:-1]
        pol=np.polyfit(t,theta,2)
   
    elif d==1:
        theta=thetar[-3:-1]
        t=tr[-3:-1]
        pol=np.polyfit(t,theta,1)
    
    #finding roots of the polynomial
    root=np.polynomial.polynomial.polyroots(pol[::-1])
    
    #return root in the range 
    for i in range(0,d):
        if t[0]<root[i]<t[d]:
            return root[i]*4


def polyin2(tr,thetar,d):
#estimates a polynomial but across the whole data set & for any d
   
    import numpy as np
    
    pol=np.polyfit(tr,thetar,d)

    
    #finding roots of the polynomial
    root=np.polynomial.polynomial.polyroots(pol[::-1])
    
    #return root in the range 
    for i in range(0,d):
        if tr[-3]<root[i]<tr[-2]:
            return root[i]*4