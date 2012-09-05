"""Module for solving Equation 6 of Canto, Raga, & Wilkin 1996 (CRW96)

This finds the radius of a stellar wind interaction bowshock in terms
of the momentum and angular momentum injected by the two winds.

CRW96 concentrated on the case of isotropic winds, but this module
will work generally with any cylindrically symmetric wind

All positions and distances are in units of D, the interstar
separation

The philosophy of this module is to be as straightforward and
transparent an implementation as possible of the equations in the
paper.  There is no attempt to be efficient.

"""
import numpy as np
import scipy.integrate
import scipy.optimize

###
### Public functions 
###
def isotropic_momentum(theta):
    """Momentum as a function of angle for an isotropic wind"""
    return 1.0

DIFFUSE_BETA = 0.0              # Parameter giving relative strength of diffuse field
def proplyd_momentum(theta): 
    """Momentum as a function of angle for a proplyd wind

    Proportional to sqrt(cos(theta)) in the head (theta < pi/2), and
    then a small constant value (zero by default) in the tail (theta >
    pi/2).  The tail value is set via the module-level variable
    DIFFUSE_BETA.

    """
    return DIFFUSE_BETA + (1.0 - DIFFUSE_BETA)*np.sqrt(max(0.0,np.cos(theta)))
    


###
### Public classes
###

class Wind(object):
    """Class to represent a stellar wind (or proplyd, etc)

    axial_momentum_flux is the momentum flux at theta=0 (or mu=1) in
    arbitrary units

    momentum_law is a function that describes how the momentum flux
    per unit solid angle varies with theta

    origin is a flag that is True if the center of this wind is at the
    co-ordinate origin.  If origin is False, then the center of the
    wind is at unit distance from the origin along the z axis.

    """
    def __init__(self, axial_momentum_flux=1.0, momentum_law=isotropic_momentum, origin=True):
        self.axial_momentum_flux = axial_momentum_flux
        self.momentum_law = momentum_law
        self.origin = origin


    def Jdot(self, theta):
        """
        Angular momentum injection rate about the origin, integrated between axis and theta
        """
        if self.origin:
            return 0.0
        else:
            if self.momentum_law == isotropic_momentum:
                return 0.25*self.axial_momentum_flux*(theta - np.sin(theta)*np.cos(theta))
            else:
                # I haven't implemented the numerical integration yet
                # in this case, but hopefully we will not need it
                raise NotImplementedError


    def Pidot_z(self, theta):
        """Linear z-momentum injection rate, integrated between axis
        and theta

        """
        if self.momentum_law == isotropic_momentum:
            # Analytic solution for isotropic case
            Pdz = 0.25*np.sin(theta)**2
        else:
            # Numerical integration for the general case
            Pdz, err = scipy.integrate.quad(self._integrand_Pdz, 0.0, theta)
        if self.origin:
            return Pdz*self.axial_momentum_flux
        else:
            # The second star has oppositely directed axial momentum
            return -Pdz*self.axial_momentum_flux

    def Pidot_r(self, theta):
        """Linear r-momentum injection rate, integrated between axis
        and theta

        """
        if self.momentum_law == isotropic_momentum:
            # Analytic solution for isotropic case
            Pdr = 0.25*(theta - np.sin(theta)*np.cos(theta))
        else:
            # Numerical integration for the general case
            Pdr, err = scipy.integrate.quad(self._integrand_Pdr, 0.0, theta)
        return Pdr*self.axial_momentum_flux

    def _integrand_Pdz(self, t):
        return 0.5*np.cos(t)*self.momentum_law(t)*np.sin(t)
    def _integrand_Pdr(self, t):
        return 0.5*np.sin(t)*self.momentum_law(t)*np.sin(t)

                            
###
### Public functions 
###
def isotropic_momentum(theta):
    """
    Momentum as a function of angle for an isotropic wind
    """
    return 1.0

DIFFUSE_BETA = 0.0              # Parameter giving relative strength of diffuse field
def proplyd_momentum(theta): 
    """
    Momentum as a function of angle for a proplyd wind
    """
    return DIFFUSE_BETA + (1.0 - DIFFUSE_BETA)*np.sqrt(max(0.0,np.cos(theta)))
    

###
### Private functions
###

def _radius_eq6(w, w1, th, th1):
    """
    Literal implementation of CRW96 Eq 6 for two winds w, w1

    Returns the radius for a given pair of angles th and th1 in terms
    of the momentum rates injected by the two winds
    """
    numerator = w.Jdot(th) + w1.Jdot(th1)
    denominator = (w.Pidot_z(th) + w1.Pidot_z(th1))*np.cos(th) \
                  - (w.Pidot_r(th) + w1.Pidot_r(th1))*np.sin(th)
    return numerator/denominator


def _radius_eq23(th, th1):
    """
    Literal implementation of CRW Eq 23 

    Gives the radius in terms of the two angles th and th1
    """
    return np.sin(th1)/np.sin(th+th1)

if __name__ = "__main__":
    
    wind = Wind(momentum_law=proplyd_momentum)
    wind1 = Wind(momentum_law=isotropic_momentum, origin=False)

