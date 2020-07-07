# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 13:57:32 2017

@author: kmiri_000
"""
import subroutines as su, numpy as np, matplotlib.pyplot as plt, scipy as sc

#11 values for dt multiplied by factor of 2
dt=np.array([0.004,0.008,0.016,0.032,0.064,0.128])

#initializing arrays for the error plots
rk1=np.array([])
rk2=np.array([])
rk3=np.array([])
rk4=np.array([])

#run for all values of dt
for i in range(0,6):
    rk1=np.append(rk1,su.test_rk(dt[i],1))
    rk2=np.append(rk2,su.test_rk(dt[i],2))
    rk3=np.append(rk3,su.test_rk(dt[i],3))
    rk4=np.append(rk4,su.test_rk(dt[i],4))

plt.plot(np.log2(dt),np.log2(rk1),'bx')
plt.plot(np.log2(dt),np.log2(rk1),label="Euler")
plt.plot(np.log2(dt),np.log2(rk2),'yx')
plt.plot(np.log2(dt),np.log2(rk2),label="RK2")
plt.plot(np.log2(dt),np.log2(rk3),'gx')
plt.plot(np.log2(dt),np.log2(rk3),label="RK3")
plt.plot(np.log2(dt),np.log2(rk4),'rx')
plt.plot(np.log2(dt),np.log2(rk4),label="RK4")
plt.legend()
plt.xlabel("Value of " r'$\log_2 \Delta t$')
plt.ylabel("Value of " r'$\log_2E$')
plt.show()

print (np.log2(rk1[5])-np.log2(rk1[0]))/(np.log2(dt[5])-np.log2(dt[0]))
print (np.log2(rk2[5])-np.log2(rk2[0]))/(np.log2(dt[5])-np.log2(dt[0]))
print (np.log2(rk3[5])-np.log2(rk3[0]))/(np.log2(dt[5])-np.log2(dt[0]))
print (np.log2(rk4[5])-np.log2(rk4[0]))/(np.log2(dt[5])-np.log2(dt[0]))

print (sc.special.gamma(1.0/4))**2/np.sqrt(np.pi)