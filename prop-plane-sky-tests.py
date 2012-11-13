import numpy as np
import matplotlib.pyplot as plt
from equation6 import Shell
import argparse
from scipy.optimize import fsolve,bisect
"""
In this script, we'll create the proyection of proplyd in the sky plane, assuming the plane of proplyd is rotated respect plane of sky by an angle i
"""

def theta_lim(beta): 
    "Asymptotic opening angle from CRW eq (28)"
    def f(theta):
	"Function to be zeroed: f(theta) = 0 for theta = theta_lim"
	return theta - np.tan(theta) - np.pi/(1.0 - beta)
    return bisect(f, 0.5*np.pi+0.01, np.pi)

def OMEGA(r,rpdr,h):

    dR = (rpdr - r)/h
    w = dR/r
    return w

def omega(r,t):
    """
    Will's version of omega(theta), looks like works better
    uses d(log(r))/dtheta = 1/r * dr/dtheta (ingenious)
    """
    womega = np.zeros_like(r)
    womega[:-1] = np.diff(np.log(r))/np.diff(theta)
    return womega

Nth = 200
Ninc = 10

beta = np.array([0.16, 0.08, 0.04, 0.02, 0.01])
innertype = ['isotropic','proplyd']
lw = dict(isotropic = 1, proplyd = 2)
for b in beta:
    print "******beta = {}******".format(b)
    thlim = theta_lim(b) 
    theta = np.linspace(0,thlim,Nth)
    for inn in innertype:
        if inn == 'proplyd':
            inc = np.linspace(0.0, 0.5*np.pi, Ninc)   
        else:
            inc = np.linspace(0.,0.98*(thlim - 0.5*np.pi),Ninc)
        print "Maximum Inclination: ", np.degrees(inc[-1])
        shell = Shell(beta=b,innertype=inn)
#       Calculing R(theta)
        R=shell.radius(theta)
        R[R<=0] = np.nan # Set neagtive R to NaN
        w = omega(R,theta)
        for jnc in inc:
            SenodePhiT = np.tan(jnc) * ( ( 1+ w*np.tan(theta) )/( w-np.tan(theta) ) )
            SenodePhiT[np.abs(SenodePhiT)>=1.] =np.nan 
            xi = R*np.cos(theta)*np.cos(jnc)-R*np.sin(theta)*SenodePhiT*np.sin(jnc)
            yi = R*np.sin(theta)*np.sqrt(1-SenodePhiT**2) 
            if inn == "isotropic":
                label = 'i={}'.format(jnc)
            else:
                label = None
            plt.plot(xi,yi,linewidth=lw[inn],label=label)                             
#    plt.axis([-1.0,1.0,-0.05,yi[-1]+0.5])
    plt.legend()
    plt.xlabel('z')
    plt.ylabel('r')
    plt.title('Bowshock shapes for wind and beta = {}'.format(b))
    plt.savefig('bowshock-wind-plane-sky-beta-{}.pdf'.format(b))
    plt.clf()
