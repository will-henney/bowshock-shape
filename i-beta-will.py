"""
La idea es crear un script que tenga una funcion que me de Rc'\D' y R0'/D' 
dandole como entrada beta, para toda i posible. Luego, dar un intervalo para 
Rc/Ro y para R0/D y de la salida que me dio ver cuales valores de i dan un R0 
y Rc que caigan en el intervalo, y guardar esos valores de i en un archivo, y 
graficar despues  
"""
#import libprop
import numpy as np
import argparse
import matplotlib.pyplot as plt
from scipy.optimize import bisect, leastsq
from equation6 import Shell
import json
from scipy.interpolate import interp1d
"""
* First step:  Set interval for R0 and Rc
* Second step: Create the bowshock NA
* Third step:  Measure Ro and Rc NA
* 2nd A step: Read the file with de R0-Rc curves
Repeat second A -and third step- for all possible inclinations (End of first test)
* Fourth step: Check if R0 and Rc fits into the interval (Possible results: Yes or No)
* Fifth step: if forurth is true, graph the corresponding b & i, if false, don't do anything
* Sixth step: Repeat from second to fifth step for many i's and b's
"""

"""
set interval:
For LV3, for example:
R0/D = 0.336 +- 1e-3
Rc/R0 = 2.0342262 +- 0.139
And so on for the rest of the proplyds
"""


class Proplyd(object):
    def __init__(self, name, beta=(0.01, 0.001), inc=(30.0, 15.0), color="r"):
        self.name = name
        self.beta, self.dbeta = beta
        self.inc, self.dinc = inc
        self.color = color



shelldata = json.load(open("rc-r0.json"))


proplyds = [
    # inclination from Henney et al 2002
    Proplyd("LV2", beta=(0.126, 0.01), inc=(40.0, 10.0), color="r"),
    Proplyd("LV3", beta=(0.061, 0.015), inc=(45.0, 15.0), color="g"),
    Proplyd("LV4", beta=(0.126, 0.01), inc=(40.0, 10.0), color="b"),
    Proplyd("LV5", beta=(0.126, 0.01), inc=(40.0, 10.0), color="y"),
    Proplyd("177-341", beta=(0.126, 0.01), inc=(40.0, 10.0), color="c"),
    Proplyd("167-328", beta=(0.126, 0.01), inc=(40.0, 10.0), color="m"),
]

#input of observational measurements of R0/D
proplyd = ["LV2","LV3","LV4","LV5","177-341","167-328"]
color = ['r','g','b','y','c','m']
obs_beta = [0.126,0.061,0.040,0.073,0.135, None]
del_beta = [0.010,0.015,0.007,0.014,0.021,None]
obs_inc = [60,45,45,45,60, None]
mirror_inc = [30,45,45,45,30,None] #in the GAH 2002 data the reported inclination is the complementary angle of
                          #the inclinations in my model
del_inc = [7,15,15,15,7, None]

# Will's original changes - now superseded
# obs_beta = [0.126,0.061,0.040,0.073,0.135, None] 
# obs_inc = [40,45,45,45,10, None]
# del_inc = [10,15,15,15,5, None]

R0m   = np.array([0.2385,0.336,0.188,0.2125,0.132,0.096])

#input of observational inputs for Rc and errorbars
Rcm   = np.array([1.468,2.034,1.987,1.501,1.405,1.297])
Rcd   = np.array([0.194,0.139,0.072,0.146,0.118,0.269])
ylow  = Rcm-Rcd  
yhigh = Rcm+Rcd

R0_epsilon = 0.1               # Assume 10% error in radii


for j,p in enumerate(proplyd):
    print p
    label = p
    # Plot the beta-inc derived from HA98 parameters
    if obs_beta[j] is not None:
        plt.errorbar(obs_inc[j], obs_beta[j], xerr=del_inc[j],yerr = del_beta[j], fmt=color[j]+'o', label="")
        plt.errorbar(mirror_inc[j], obs_beta[j], xerr=del_inc[j],yerr = del_beta[j], fmt=color[j]+'D', label="")
        #Also plot for the complementary inclinations
    for beta, beta_data in shelldata.items():
        beta = float(beta)
        r0 = np.array(beta_data["R0'"])
        rc = np.array(beta_data["Rc"])/r0
        inc = np.array(beta_data['inc'])

        # Select all points consistent with error bars
        m1 = abs(rc - Rcm[j]) <= Rcd[j] # Rc'/R0' axis
        m2 = abs(r0 - R0m[j]) <= R0_epsilon*R0m[j] # R0'/D' axis
        inc_good = inc[m1 & m2] # Points must satisfy both conditions
        ngood = len(inc_good)
        beta_good = np.ones((ngood,))*beta
        if ngood > 0:
            plt.plot(np.degrees(inc_good), beta_good, color[j]+'.', label=label, alpha=1.0)
            label = ""
            
        # Also plot points that only agree with R0'/D'
        inc_good = inc[m2] 
        ngood = len(inc_good)
        beta_good = np.ones((ngood,))*beta
        if ngood > 0:
            plt.plot(np.degrees(inc_good), beta_good, color[j]+'.', label="", alpha=0.1)


plt.yscale('log')            
plt.grid()
plt.xlim(0,90)
plt.ylim(0.001 - 1e-4, 0.16 + 1e-4)
plt.xlabel("i(deg)")
plt.ylabel("beta")
plt.legend(loc="best")
plt.title("i vs beta plausible for proplyds")
plt.savefig("i-beta-will.pdf")
