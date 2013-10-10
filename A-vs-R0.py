"""
The goal of this program is to replicate the final graph made in protocolo 
using the analytic functions found by Will for A projected in function of R0/D 
projected
"""

import numpy as np
from matplotlib import pyplot as plt

#1
# Define my beta values and inclination angles

beta = np.array([0.001,0.005,0.01,0.05,0.08])
betacolor = ["r","g","b","c","m"]
inc = np.radians([0.0,15.0,30.0,45.0,60.0,75.0])
#2
# Compute R0/D and A as a function of beta
innertype = ["isotropic","proplyd"]

y0 = {"proplyd":0.66,"isotropic": 0.585}
b = {"proplyd":50.0,"isotropic":5.0}
d = {"proplyd": 0.49, "isotropic": 0.5}
opacity = {"isotropic":0.2,"proplyd":1.0}

for inn in innertype:
    c = np.exp(-1./b[inn])
    for i,B in enumerate(beta):
        R0 = np.sqrt(B)/(1+np.sqrt(B))
        A =  (1-c)/(y0[inn]*(np.exp(-B**d[inn]/b[inn])-c))
        a = (A-1)/A
#3 
# Use Will notes to determine the projected A and R0 for different inclination angles

        R0p = R0*A*(1-a*np.cos(inc))/np.cos(inc)
        Ap = 1./(1-a*np.cos(inc))

#4 
# Plot plot plot
        label = {"isotropic":None,"proplyd":"beta={}".format(B)}
        plt.plot(R0p,Ap,color=betacolor[i],alpha=opacity[inn],label=label[inn])

plt.legend(loc="best")
plt.grid()
plt.xlabel("R0'/D'")
plt.ylabel("Rc'/R0'")
#plt.axis("equal")
plt.xlim(0.0,0.5)
plt.ylim(0,4.)
plt.savefig("A-vsR0.pdf")
