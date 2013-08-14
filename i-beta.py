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
"""

R0m   = 0.336
R0d   = 1e-3
xlow  = R0m-R0d
xhigh = R0m+R0d

Rcm   = 2.0342262
Rcd   = 0.139
ylow  = Rcm-Rcd  
yhigh = Rcm+Rcd


shelldata = json.load(open("rc-r0.json"))

for model in shelldata.items():
    uni_beta,radii = model
    beta = float(uni_beta)
    r0 = np.array(radii["R0'"])
    rc = np.array(radii['Rc'])
    inc = np.array(radii['inc'])
    
# choose the matching radii with observations
    m1 = (r0 < xhigh) & (r0 > xlow)
    m2 = (rc/r0 < yhigh) & (rc/r0 > ylow)
# if the m1 and m2 conditions are satisfied, then plot the data, the y axis is beta and the
# x axis is the inclination. This only applies for LV3 so far
    plt.plot(np.degrees(inc[m1 & m2]),np.ones(len(inc[m1&m2]))*beta,'r.')

plt.xlim(0,90)
plt.ylim(0.001 - 1e-4,0.08 + 1e-4)
plt.xlabel("i(deg)")
plt.ylabel("beta")
plt.title("i vs beta plausible for LV3")
plt.savefig("i-beta.pdf")
