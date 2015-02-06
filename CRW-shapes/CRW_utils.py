import numpy as np
from scipy.interpolate import interp1d
class projection(object):
"""
Class to obtain the projected shape of a given 
shape with cilindrical symetry
"""
    def __init__(self,r,theta,i):

        self.r = r
        self.theta = theta
        self.i = np.radians(i)
    
    def omega(self):
        w = np.zeros_like(self.r)
        w[:-1] = np.diff(np.log(self.r))/np.diff(self.theta)
        w[-1] = w[-2]
        return w

    def tangent_alpha(self):

        w = self.omega()
        return (1+w*np.tan(self.theta))/(np.tan(self.theta)-w)

    def tangent_phi(self):
        ta = self.tangent_alpha()
        tp = -ta*np.tan(self.i)
        if np.abs(tp) > 1.0:
            return np.nan
        else:
            return tp


    def rotation(self):
        tphi = self.tangent_phi()
        xpt = self.r*(np.cos(self.theta)-np.sin(self.theta)*tphi*np.tan(self.i))
        ypt = self.r*np.tan(self.theta)*np.sqrt(1-tphi**2)
        zpt = self.r*(np.cos(self.theta)*np.sin(self.i)+np.sin(theta)*tphi*np.cos(self.i))
        return xpt,ypt,zpt


    def spher_primed(self):
        xp,yp,zp = self.rotation()
        rp = np.sqrt(xp**2+yp**2+zp**2)
        tp = np.arctan2(yp,xpq)
        return rp,tp


class char_radii(object):
"""
Class to obtain characteristic radii (Regardless projection)
Input: R',theta'
"""
    def __init__(self,R,t):

        self.R = R
        self.t = t

    def R0(self):
        return R.min()

    def R90(self):

        f = interp1d(self.t,self.R,bounds_error=False)
        return f(0.5*np.pi)
        

    def Rc(self):
        """
        I only need the radius of curvature at theta'=0,
        so I will calculate the derivatives at that point
        """
        rdot = np.diff(self.R)/np.diff(self.theta)
        rddot = np.diff(rdot)/np.diff(self.theta[:-1])

        return (self.R0()**2+rdot[0]**2)**(1.5)/np.abs(self.R0()**2+2*rdot[0]**2-self.R0()*rddot[0])


class thickness(object):
    """
    Computes the thickness of the bowshock
    """

    def __init__(self,M,r,t):
        
        self.M = M
        self.r = r
        self.t = t

    def inner_r(self):
        
        
    
