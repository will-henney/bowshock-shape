"""
In this script, we create the projection of proplyd in the sky
plane, assuming the proplyd axis is at an angle i with respect to
plane of sky
"""

import argparse
import json
import conic_utils
import numpy as np
from scipy.optimize import bisect
from scipy.interpolate import interp1d
import lmfit
from equation6 import Shell


#
# Command line arguments
#
parser = argparse.ArgumentParser(
    description="""Find characteristic radii for projected 
    bowshock models""")

parser.add_argument("--mach","-m",type=float,default=3.0,help="mach number of inner wind")
cmd_args = parser.parse_args()


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
    womega[:-1] = np.diff(np.log(r))/np.diff(theta)
    return womega





#First, design a set of beta values to calculate R(theta) (Well, first
#than that, desgign a theta array)


# beta = [0.08, 0.04, 0.02, 0.01, 0.005, 0.002, 0.001]
#beta = [0.005, 0.002, 0.001]
beta = np.linspace(0.08,0.001,50)
Nth = 800
Ninc = 100

innertype = ["isotropic", "proplyd"]
# innertype = ['isotropic']

shelldata = {"isotropic":"","proplyd":""}
#Given the mach number, we can calculate the index of \cos\psi in the shell thickness and
#the thickness at the axis
H0 = 3./(4*cmd_args.mach**2+1.)
n  = 3.4*(4.*cmd_args.mach**2+1)/(4.*cmd_args.mach**2+35) 
for inn in innertype:
    shelldat={}
    print "******{} case******".format(inn)
    for b in beta:
        print "beta = ", b
        thlim = theta_lim(b)  # I have to calculate theta_lim for
                             # proplyd case
        theta = np.linspace(0, thlim, Nth)
        if inn == 'proplyd':
            # In proplyd case, try all inclinations
            inclinations = np.linspace(0.0, 0.5*np.pi, Ninc)
        else:
            inclinations = np.linspace(0., 0.98*(thlim - 0.5*np.pi), Ninc)
            print "Maximum Inclination: ", np.degrees(inclinations[-1])

        # Initialize the data arrays to NaNs. R0 and R90 are enough to calculate
        # the a,b parameters of an ellipse or hyperbola
        R0 = np.zeros_like(inclinations)*np.nan
        R90 = np.zeros_like(inclinations)*np.nan
        thinf = np.zeros_like(inclinations)*np.nan
        
        shell = Shell(beta=b, innertype=inn)
        R = shell.radius(theta)
        w = omega(R, theta)
        tan_alpha = (1 + w*np.tan(theta)) /   (w - np.tan(theta))
        alpha = -np.arctan(tan_alpha)
        psi = theta+alpha-0.5*np.pi
        psi_l = np.arccos(1./cmd_args.mach)
        psi[psi>psi_l] = np.nan
        H = H0*np.cos(psi)**(-n)
        Rin = R- H*R[0]/np.cos(psi)
        w_in = omega(Rin,theta)
        tan_alpha_in = (1+w_in*np.tan(theta))/(w_in-np.tan(theta))
        for i, inc in enumerate(inclinations):
            SenodePhiT = np.tan(inc)*tan_alpha_in 
            # other way to set mask
            SenodePhiT[np.abs(SenodePhiT) >= 1.] = np.nan

            # Correct for projection and for fact that observed radii are normalised by
            # D' = D cos(inc) (Now for inner radius)
            xi = (Rin/np.cos(inc))*(np.cos(theta)*np.cos(inc)
                                  - np.sin(theta)*SenodePhiT*np.sin(inc))
            yi = (Rin/np.cos(inc))*np.sin(theta)*np.sqrt(1-SenodePhiT**2)
            mask = np.isfinite(xi) & np.isfinite(yi)
            xim, yim = xi[mask], yi[mask]  # Removing nan elements from xi
                                           # and yi
            Rp = np.hypot(xim,yim)
            thp = np.arctan2(yim,xim)
            #Characterictic radii of shell
            try:
                R0[i] = xim[0]
            except IndexError:
                print("max inclination={}".format(inc))
                break
            try:
                f=interp1d(thp,Rp,bounds_error=False) 
                R90[i] = f(0.5*np.pi)
            except ValueError:
                print("max inclination={}".format(np.degrees(inc)))   
                break   
            #Conic parameters got from R0 and R90
            if 2*R0[i] > R90[i]:
                #Ellipse case
                a = R0[i]**2/(2*R0[i]-R90[i])
                b2 = np.sqrt(a*R90[i])
                thinf[i] = np.degrees(np.arctan2(b2,a))
            elif 2*R0[i] == R90[i]:
                #parabola case
                thinf[i] == 0.0
            else:
                # hyperbola case
                a = np.abs(R0[i]**2/(R90[i]-2*R0[i])) 
                b2 = np.sqrt(a*R90[i])
                thinf[i] = -np.degrees(np.arctan2(b2,a))
        shelldat[b]={
          "inc": inclinations[np.isfinite(R0) & np.isfinite(R90)].astype(float).tolist(),
          "R0": R0[np.isfinite(R0) & np.isfinite(R90)].astype(float).tolist(),
          "R90": R90[np.isfinite(R0) & np.isfinite(R90)].astype(float).tolist(),
          "thinf": thinf[np.isfinite(R0) & np.isfinite(R90)].astype(float).tolist()
        }
    shelldata[inn]=shelldat
# Save all the theoretical curves
with open("r90-r0-in-{}.json".format(cmd_args.mach), "w") as f:
    json.dump(shelldata,f,indent=2)
