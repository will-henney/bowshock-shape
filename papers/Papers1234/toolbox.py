"""
Set of functions to calculate
tangent line, projected shapes, 
characteristic radii, etc.
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import leastsq
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

def R0p(R,t):
    """
    Returns R'_0 and a uncertainty due to the first element of theta' is not
    (almost sure) exactly 0 
    """


    r0p = R[0]
    w = omega(R,t)
    dr = np.abs(r0p*w[0]*t[0])
    return r0p,dr

def Rc(R,t):
    """
    Computes the radius of curvature at the shell's <<nose>> and the error margin.
    This is a numeric approach, so, for i=0, the result will differ from the analytic 
    result, but is accurate enough. Note: Enter projected R & theta
    """
    w = omega(R,t)
    dw = np.gradient(w)/np.gradient(t) # Derivative of omega
    rc = R[0]/np.abs(-dw[0]+1)
    f = (w**2+1)**1.5/np.abs(w**2-dw+1)
    df = np.gradient(f)/np.gradient(t)
    drc = np.abs(R[0]*(df[0])*t[0])
    return rc,drc

def Rc_fit(R,t):
    """
    Compute Rc via least squares fit
    (the old fashion way)
    """
    def complete_curve(R,t):
        """
        return the complete shape of the bowshock in cartesian coordinates
        """
        theta_com = np.concatenate((sorted(-t),t))
        R_com = np.concatenate((sorted(R,reverse=True),R))
        x_com,y_com = R_com*np.cos(theta_com),R_com*np.sin(theta_com)
        return x_com,y_com
    
    def R_calc(xc,yc):
        """
        calculate the distance between the (x,y) curve and the center (xc,yc)
        """
        return np.sqrt((xfit-xc)**2+(yfit-yc)**2)
    def f_fit(c):
        Ri = R_calc(*c)
        return Ri-Ri.mean()
    x,y = complete_curve(R,t)
    m = np.abs(np.degrees(np.arctan2(y,x)))-45. <=0
    xfit,yfit = x[m],y[m]
    c0 = 0,0
    c1,ier = leastsq(f_fit,c0)
    xf,yf = c1
    ri = R_calc(xf,yf)
    return ri.mean()

def R90(R,t):
    """
    Computes the radius of the shell perpendicular to the symmetry axis. And the 
    corresponding errorbar. Insert the projected R and theta.
    """

    r_int = interp1d(t,R,bounds_error= False)
    R90 = r_int(0.5*np.pi)
    return R90
