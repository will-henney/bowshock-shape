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
    w[:-1]=np.diff(np.log(R))/np.diff(t)
    w[-1]=w[-2] #Assumes that the slope of the tangent line
                # in the wings is nearly constant
    return w

def tangent_line(i,w,t):
    """
    sine of the azimutal angle of the tangent line
    """
    return np.tan(np.radians(i))*(1+w*np.tan(t))/(w-np.tan(t))

def projection(R,t,i,sp):
    """
    Coordinates of the tangent line in the observer's
    reference frame
    """
    xp = R*(np.cos(t)-np.sin(t)*sp*np.tan(np.radians(i)))
    yp = R*(np.sin(t))*np.sqrt(1-sp**2)/np.cos(np.radians(i))
    return xp,yp

def theta_prime(xp,yp):
    """
    Polar angle of the projected shell
    """
    return np.arctan2(yp,xp)

def R_prime(xp,yp):
    """
    Returns projected radius of the shell
    """
    return np.sqrt(xp**2+yp**2)

def Rc(R0,w,tp):
    """
    Computes the radius of curvature at
    the shell's <<nose>>
    """
    dw = (w[1]-w[0])/tp[1]-tp[0] # Derivative of omega at theta'=0
    return r0*(w[0]**2+1)**1.5/np.abs(w[0]**2-dw+1)
