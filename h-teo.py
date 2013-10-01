from equation6 import Shell
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import bisect
"""
Main goal: Trace the theoretical inner and outer shell in bowshock
Hypothesis:
1.-  v*cos(\beta_a)/c \equiv M_\perp.
 Where v is the inner wind velocity, 
c is the sound speed, beta_a is the angle between the velocity vector and the normal
 vector in the shell and M_\perp is the mach number perpendicular to the shell.
2.- h \propto (M_\perp)**{-2}.
Where h is the shell width.
So, we conclude that h/h_0 = (cos(\beta_a))**{-2} where h_0 is the shell width at theta=0
3.- R_in = R - h
"""

#Steps:
#1: Create a theta array and create a bowshock shell using CRW formalism
#2: Calculate \beta_a for each theta: we can use geometrical arguments
#3: Calculate h and R_in
#4: Plot everything

def theta_lim(beta):
    "Asymptotic opening angle from CRW eq (28)"
    def f(theta):
        "Function to be zeroed: f(theta) = 0 for theta = theta_lim"
        return theta - np.tan(theta) - np.pi/(1.0 - beta)
    return bisect(f, 0.5*np.pi+0.01, np.pi)

def omega(r, t):
    """
    Will's version of omega(theta), looks like work better.  uses
    d(log(r))/dtheta = 1/r * dr/dtheta (ingenious)
    """
    womega = np.zeros_like(r)
    womega[:-1] = np.diff(np.log(r))/np.diff(t)
    return womega

#1:
Nth = 800
th_lim = theta_lim(0.01)
#theta = np.linspace(0,th_lim,Nth)
shell = Shell(beta=0.01,innertype="proplyd")
#R_ext = shell.radius(theta)
#w = omega(R_ext,theta) 
#for r,t,m in zip(R_ext,theta,w):
#    print r,np.degrees(t),m

#2: using geometrical arguments, we can say that 
# \beta_a = alpha + theta - 90 so
# (\cos(\beta_a))**2 = (\sin(alpha+theta))**2
# where \tan(alpha) = (1+w*\tan(theta))/(\tan(theta)-w) 

#3: switch to circular approximation for calculating R_in
# First: design an angle theta_c

#Add the analytic fit found for A = Rc/R0
y0 = 0.66
b = 50.0
c = np.exp(-1./b)
d = 0.49
A =  (1-c)/(y0*(np.exp(-0.01**d/b)-c))
a = (A-1)/A
theta_c = np.linspace(0,0.5*np.pi,Nth)
theta = np.arctan2(np.sin(theta_c),np.cos(theta_c)-a)

R_ext = shell.radius(theta)
R_circ = R_ext[theta==0]*np.sqrt( 1+( 2*a*( 1-np.cos(theta_c) ) ) / (1-a)**2 )

#alpha = np.arctan2(1+w*np.tan(theta),np.tan(theta)-w)
h0 = 0.2
sqcosine = (np.cos(theta_c-theta))**2
#h = h0*(np.sin(alpha+theta))**(-2)
h = h0/sqcosine**2
R_in = R_ext*(1-h) 

#4 Translate to cartesian coordinates and plot

x_e,y_e = R_ext*np.cos(theta),R_ext*np.sin(theta)
x_i, y_i = R_in*np.cos(theta),R_in*np.sin(theta)
x_c,y_c = R_circ*np.cos(theta),R_circ*np.sin(theta)
plt.plot(x_e,y_e,"r-",label="Outer shock (CRW formalism)")
plt.plot(x_i,y_i,"g-",label="Inner Shock")
plt.plot(x_c,y_c,"b-",label="Circular approximation")
plt.legend(loc="best")
plt.axis("equal")
#plt.xlim(-0.5,0.2)
plt.ylim(0,0.25)
plt.grid()
#plt.plot(np.degrees(theta_c),h)
plt.savefig("test.pdf")
plt.clf()
