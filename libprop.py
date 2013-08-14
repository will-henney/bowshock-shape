"""
Basically, this is a library of functions I will need
"""

import argparse
import json

import numpy as np
from scipy.optimize import bisect, leastsq

from equation6 import Shell


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


def theta_lim(beta):
    "Asymptotic opening angle from CRW eq (28)"
    def f(theta):
        "Function to be zeroed: f(theta) = 0 for theta = theta_lim"
        return theta - np.tan(theta) - np.pi/(1.0 - beta)
    return bisect(f, 0.5*np.pi+0.01, np.pi)


def shellmaker(b,inc):
    """
    Creates the projected and normalized bowshock shape
    Entries:
    b: momentum rate ratio between proplyd and massive star. Scalar
    inc: orientation angle of proplyd: scalar
    """
    Nth = 800
    tfit = 45
    
    shell = Shell(beta=b, innertype="proplyd")
    thlim = theta_lim(b)
    theta = np.linspace(0, thlim, Nth)
    R = shell.radius(theta)
    w = omega(R, theta)
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
    return x2,y2

def measure_radii(b,inc):

    def calc_R(xc, yc):
        """ calculate the distance of each 2D points from the center (xc, yc) """
        return np.sqrt((xm-xc)**2 + (ym-yc)**2)


    def f_2(c):
        """Calculate the algebraic distance between the 2D points
        and the mean circle centered at c=(xc, yc)"""
        Ri = calc_R(*c)
        return Ri - Ri.mean()


    #initialize R0 and Rc to NaN
    R0 = np.zeros_like(inc)*np.nan
    Rc = np.zeros_like(inc)*np.nan
    # initial guess for the circle fit to calculate curvature radius
    c0 = 0, 0

    for i, incl in enumerate(inc):
        xs,ys = shellmaker(b,incl)
        m = np.abs(np.degrees(np.arctan2(ys, xs))) - 45 <= 0.
        print i
	try:
            xm,ym = xs[m],ys[m]
	except IndexError:
	    print "IndexError, max inclination = ",incl
	    break
        try:
            c_fit, ier = leastsq(f_2, c0)
            xfit, yfit = c_fit
            Ri_2 = calc_R(xfit, yfit)
            R0[i] = xm[0]
        # The fit uses the bowshock data such that t<45 deg
            Rc[i] = Ri_2.mean()
        except TypeError:
            print 'error'
	    break
        #skip invalid inputs for fitting and skip to next beta value
    return R0,Rc
