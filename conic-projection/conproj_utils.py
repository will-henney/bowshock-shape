"""Projection of conic sections (ellipses ... hyperbolae)

"""
from __future__ import print_function
import numpy as np

class Conic(object):
    """Bowshock shape - surface of revolution of a plane conic section

    As the shape parameter `th_conic` is varied, this gives a sequence
    from (`th_conic` = 45 - 90) oblate spheroid -> (`th_conic` = 45)
    sphere -> (`th_conic` = 0 - 45) prolate spheroid -> (`th_conic` =
    0) paraboloid -> (`th_conic` = -90 - 0) hyperboloid.

    All the shapes are normalized so that the nose of the bow is at
    unit distance from the bowshock source (star or proplyd), along
    the x-axis.  The radius of curvature on this axis is set by
    another shape parameter, `A`

    The parameter `t` is not an angle, but the angle can be found with
    :method:`theta`

    """
    tlimits_e = -np.pi, np.pi
    tlimits_h = -5.0, 5.0
    def __init__(self, A=1.0, th_conic=45.0):
        # Hyperbolic or not?
        self.hyper = th_conic < 0.0
        # Axis ratio
        self.b_a = np.tan(np.radians(abs(th_conic)))
        # Scaled radius of curvature: Rc/r0
        self.A = A

    def make_t_array(self, limits=None, N=2001):
        """Return an array of `N` equal-spaced t-parameter values between `limits`"""
        if limits is None:
            limits = self.tlimits_h if self.hyper else self.tlimits_e
        return np.linspace(limits[0], limits[1], N)

    def x(self, t):
        """Body-frame x-position wrt star"""
        fac = 1.0 - np.cosh(t) if self.hyper else np.cos(t) - 1.0
        return 1.0 + self.A*fac/self.b_a**2

    def y(self, t):
        """Body-frame y-position wrt star"""
        fac = np.sinh(t) if self.hyper else np.sin(t)
        return self.A*fac/self.b_a

    def theta(self, t):
        """Body-frame angle from x-axis in degrees"""
        return np.degrees(np.arctan2(self.y(t), self.x(t)))

    def xt(self, inc, t):
        """Observer-frame x'-position of tangent line"""
        fac1 = 1.0 - np.cosh(t) if self.hyper else np.cos(t) - 1.0
        fac2 = np.cosh(t) if self.hyper else np.cos(t)
        cosi = np.cos(np.radians(inc))
        tani = np.tan(np.radians(inc))
        return cosi*(1.0 + self.A*(fac1/self.b_a**2 + fac2*tani**2))

    def yt(self, inc, t):
        """Observer-frame y'-position of tangent line"""
        st = np.sinh(t) if self.hyper else np.sin(t)
        ct = np.cosh(t) if self.hyper else np.cos(t)
        tani = np.tan(np.radians(inc))
        return np.sign(t)*(self.A/self.b_a)*np.sqrt(st**2 - (self.b_a*tani*ct)**2)

    def theta_t(self, inc, t):
        """Observer-frame angle of tangent line from x'-axis in degrees"""
        return np.degrees(np.arctan2(self.yt(inc, t), self.xt(inc, t)))

    def tparallel(self, inc):
        """Minimum value of t on the tangent line as function of inclination

        Corresponds to y'=0"""
        at = np.arctanh if self.hyper else np.arctan
        tani = np.tan(np.radians(inc))
        return at(self.b_a*tani) 
        
if __name__ == "__main__":
    # TODO Test the mechanism by drawing some conics
    print("Write some tests!")
