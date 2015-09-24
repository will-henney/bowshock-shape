"""
Set of functions to calculate
tangent line, projected shapes, 
characteristic radii, etc.
"""

import numpy as np

def omega(R,t):
    """
    Calculates the deivate of ln(R(theta))
    """
    w = np.zeros_like(t)
    w[:-1] = np.diff(np.log(R))/np.diff(t)
    w[-1]= w[-2] # assumes w is nearly constant at the shell's wings
    return w

def tangent_line(i,t,R):
    """
    sine of the azimutal angle of the tangent line
    """
    w = omega(R,t)
    return  np.tan(np.radians(i))*(1+w*np.tan(t))/(w-np.tan(t))

def projection(R,t,i):
    """
    Coordinates of the tangent line in the observer's
    reference frame
    """
    sp = tangent_line(i,t,R)
    xp = R*(np.cos(t)-np.sin(t)*sp*np.tan(np.radians(i)))
    yp = R*(np.sin(t))*np.sqrt(1-sp**2)/np.cos(np.radians(i))
    mask = np.isfinite(xp) & np.isfinite(yp)
    return xp[mask],yp[mask]

def theta_prime(xp,yp):
    """
    Polar angle of the projected shell. Enter x' & y'
    """
    return np.arctan2(yp,xp)

def R_prime(xp,yp):
    """
    Returns projected radius of the shell. Enter x' and y'
    """
    return np.sqrt(xp**2+yp**2)

def Rc(R,t):
    """
    Computes the radius of curvature at the shell's <<nose>>
    This is a numeric approach, so, for i=0, the result will differ from the analytic 
    result, but is accurate enough. Note: Enter projected R & theta
    """
    w = omega(R,t)
    dw = (w[1]-w[0])/(t[1]-t[0]) # Derivative of omega at theta'=0
    return R[0]*(w[0]**2+1)**1.5/np.abs(w[0]**2-dw+1)
