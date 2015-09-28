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

    return np.gradient(np.log(R))/np.gradient(t)

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
    dw = np.gradient(w)/np.gradient(t) # Derivative of omega
    return R[0]*(w[0]**2+1)**1.5/np.abs(w[0]**2-dw[0]+1)
