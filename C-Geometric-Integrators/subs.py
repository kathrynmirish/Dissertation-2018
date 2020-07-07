# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 11:31:44 2018

@author: kmiri_000
"""
#euler time step code
def euler(h,r,v):
    import numpy as np
    
    #defining empty array for new r
    rplus1=np.array([0.0,0.0])
    
    #defining empty array for new v
    vplus1=np.array([0.0,0.0])
    
    #updating the values for new r
    rplus1[0]=r[0]+h*v[0]
    rplus1[1]=r[1]+h*v[1]
    
    #defining unit vector rhat and constant rconst
    rconst=np.linalg.norm(r)
    rhat=r/rconst
    
    #updating the values for new v
    vplus1[0]=v[0]-(h*rhat[0])/(rconst**2)
    vplus1[1]=v[1]-(h*rhat[1])/(rconst**2)
    
    #return one array thats r and v 
    return np.append(rplus1,vplus1)

#returns analytic solution
def analytical(e,E):
    import numpy as np
    
    #define b 
    b=np.sqrt(1-e**2)
    
    #analytical solution
    r=np.array([np.cos(E)+e,b*np.sin(E)])
    
    return r

#modified euler time step code   
def modeuler(h,r,v):
    import numpy as np
    
    #defining empty array for new r
    rplus1=np.array([0.0,0.0])
    
    #defining empty array for new v
    vplus1=np.array([0.0,0.0])
    
    #updating the values for new r
    rplus1[0]=r[0]+h*v[0]
    rplus1[1]=r[1]+h*v[1]
    
    #defining unit vector rhat and constant rconst
    rconst=np.linalg.norm(rplus1)
    rhat=rplus1/rconst
    
    #updating the values for new v
    vplus1[0]=v[0]-(h*rhat[0])/(rconst**2)
    vplus1[1]=v[1]-(h*rhat[1])/(rconst**2)
    
    #return one array thats r and v 
    return np.append(rplus1,vplus1)

#leapfrog timestep code    
def leapfrog(h,r,v):
    import numpy as np
    
    #defining empy array for r'
    rprime=np.array([0.0,0.0])
    
    #defining empty array for new r
    rplus1=np.array([0.0,0.0])
    
    #defining empty array for new v
    vplus1=np.array([0.0,0.0])
    
    #updating values for r'
    rprime[0]=r[0]+(h/2)*v[0]
    rprime[1]=r[1]+(h/2)*v[1]

    #defining unit vector r'hat and constant r'
    rconst=np.linalg.norm(rprime)
    rhat=rprime/rconst    
    
    #updating values for v
    vplus1[0]=v[0]-(h*rhat[0])/(rconst**2)
    vplus1[1]=v[1]-(h*rhat[1])/(rconst**2)
    
    #updating values for r
    rplus1[0]=rprime[0]+(h/2)*vplus1[0]
    rplus1[1]=rprime[1]+(h/2)*vplus1[1]
    
    #return one array thats r and v 
    return np.append(rplus1,vplus1)

def RK4(h,r,v):
    import numpy as np

    #defining the constants
    k1=h*v
    k11=h*(-1/(np.linalg.norm(r)**3))*r
    
    k2=h*(v+0.5*k11)
    k22=h*(-1/(np.linalg.norm(r+0.5*k1)**3))*(r+0.5*k1)
    
    k3=h*(v+0.5*k22)
    k33=h*(-1/(np.linalg.norm(r+0.5*k2)**3))*(r+0.5*k2)
    
    k4=h*(v+k33)
    k44=h*(-1/(np.linalg.norm(r+k3)**3))*(r+k3)
    
    #defining new r value
    rplus1=r+((k1+2*k2+2*k3+k4)/6.0)
    
    #defining new v value
    vplus1=v+((k11+2.0*k22+2.0*k33+k44)/6.0)
    
    #return one array thats r and v 
    return np.append(rplus1,vplus1)


    
    