"""
In this script, we create the projection of proplyd in the sky
plane, assuming the proplyd axis is at an angle i with respect to
plane of sky
"""

import argparse
import json

import numpy as np
from scipy.optimize import bisect, leastsq

from equation6 import Shell


#
# Command line arguments
#
parser = argparse.ArgumentParser(
    description="""Find parallel and perpendicular projected radii for
    bowshock models""")
parser.add_argument(
    '--tobs', type=float, default=90,
    help='projected bowshock radius at theta = tobs respect to axis')
parser.add_argument(
    '--tfit', type=float, default=45,
    help='upper angle of fit data')

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


def Rt_finder(x, y, t):
    """
    This function finds a plausible value for R90, without calculating
    Theta_perp. In the rotated frame, r90 is y such that x = 0. So i
    find this value looking for the place where x has a sign change,
    if this function does not find any sign change, returns a NaN. If
    x and y don't have the same size, returns error (for completeness)
    """
    size = np.size(x)
    sizey = np.size(y)
    assert size == sizey, "x and y array sizes must be equal"

    # # If inclination is < 90, then x decreases monotonically
    # npos = len(x[x > 0.0])         # number of positive x values
    # if npos > 0 and npos > size:
    #     out = 0.5*(y[npos-1]+y[npos])
    #     return out

    for i in range(0, size-1):
        s = (np.degrees(np.arctan2(y[i], x[i])) - t)*(
            np.degrees(np.arctan2(y[i+1], x[i+1])) - t)
        if s < 0:
            slope = (y[i+1]-y[i])/(x[i+1]-x[i])
            xout = (y[i]-slope*x[i])/(np.tan(np.radians(t)) - slope)
            yout = y[i] + slope*(xout-x[i])
            return np.sqrt(xout**2+yout**2)
    return np.nan


def calc_R(xc, yc):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((xim_fit-xc)**2 + (yim_fit-yc)**2)


#@countcalls
def f_2(c):
    """Calculate the algebraic distance between the 2D points
    and the mean circle centered at c=(xc, yc)"""
    Ri = calc_R(*c)
    return Ri - Ri.mean()


#First, design a set of beta values to calculate R(theta) (Well, first
#than that, desgign a theta array)


# beta = [0.08, 0.04, 0.02, 0.01, 0.005, 0.002, 0.001]
#beta = [0.005, 0.002, 0.001]
beta = np.linspace(0.08,0.001,50)
Nth = 800
Ninc = 100

innertype = ['isotropic', 'proplyd']
# innertype = ['isotropic']

tfit = cmd_args.tfit

# Save all the curves in a dictionary for dumping
try:
    shelldata = json.load(open("rc-r0.json"))
except IOError:
    shelldata = {}

for inn in innertype:
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

        # Initialize the data arrays to NaNs
        R0 = np.zeros_like(inclinations)*np.nan
        Rt = np.zeros_like(inclinations)*np.nan
        Rc = np.zeros_like(inclinations)*np.nan

        shell = Shell(beta=b, innertype=inn)
        R = shell.radius(theta)
        w = omega(R, theta)
        for i, inc in enumerate(inclinations):
            SenodePhiT = (np.tan(inc) * ((1 + w*np.tan(theta)) /
                                         (w - np.tan(theta))))
            # other way to set mask
            SenodePhiT[np.abs(SenodePhiT) >= 1.] = np.nan

            # Correct for projection and for fact that observed radii are normalised by
            # D' = D cos(inc)
            xi = (R/np.cos(inc))*(np.cos(theta)*np.cos(inc)
                                  - np.sin(theta)*SenodePhiT*np.sin(inc))
            yi = (R/np.cos(inc))*np.sin(theta)*np.sqrt(1-SenodePhiT**2)
            mask = np.isfinite(xi) & np.isfinite(yi)
            xim, yim = xi[mask], yi[mask]  # Removing nan elements from xi
                                           # and yi
            #Creating the other part of bowshock
            xplot, yplot = np.array([xim, xim]), np.array([-yim, yim])
            x2, y2 = xplot.reshape(2*len(xim),), yplot.reshape(2*len(yim),)

            # initial guess for the circle fit to calculate curvature radius
            c0 = 0, 0
            try:
                m = np.abs(np.degrees(np.arctan2(y2, x2))) - tfit <= 0.
                xim_fit = x2[m]
                yim_fit = y2[m]
            except IndexError:
                print "Maximum Inclination:", np.degrees(inc)
                break
                #skip invalid inputs for fitting and skip to next beta value
            try:
                c_fit, ier = leastsq(f_2, c0)
                xfit, yfit = c_fit
                Ri_2 = calc_R(xfit, yfit)
                R0[i] = xim[0]
                # The fit uses the bowshock data such that t<45 deg
                Rc[i] = Ri_2.mean()
            except TypeError:
                print "Maximum Inclination:", np.degrees(inc)
                break
                #skip invalid inputs for fitting and skip to next beta value

#organize all the data in order to save everything in a file, I only need the proplyd data so:
        if inn == "proplyd":
            shelldata[b] = {
                "inc": inclinations[np.isfinite(R0)].astype(float).tolist(),
                "R0'": R0[np.isfinite(R0)].astype(float).tolist(),
                "Rc": Rc[np.isfinite(R0)].astype(float).tolist()
                }


# Save all the theoretical curves
with open("rc-r0.json", "w") as f:
    json.dump(shelldata, f)
