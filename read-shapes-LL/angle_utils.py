"""
Functions to make it easier to use astropy.coordinates objects
"""
from __future__ import print_function
from astropy import coordinates as coord
from astropy import units as u
import numpy as np

def vector_separation_old(c1, c2, mode="xy"):
    """Find vector separation between two coordinates c1, c2

    Old version, which is obsolete now that astropy.coordinates has
    more functionality

    """
    x = (c2.ra - c1.ra).to(u.arcsec).value * np.cos(c1.dec.to(u.rad).value)
    y = (c2.dec - c1.dec).to(u.arcsec).value
    if mode == "xy":
        return x, y
    elif mode == "polar":
        R = np.hypot(x, y)
        PA = np.degrees(np.arctan2(x, y) % (2*np.pi))
        return R, PA
    else:
        raise NotImplementedError

def vector_separation(c1, c2, mode="xy"):
    """Find vector separation between two coordinates c1, c2

    Input arguments c1, c2 should be astropy.coordinates objects (or
    similar).  When mode="xy" (default), return the RA and Dec
    separations in arcsec.  When mode="polar", return the scalar
    separation in arcsec and PA of the separation.  In either case,
    pure numbers are returned.  The sense of the separation is that it
    takes you from c1 to c2.

    """
    R = c1.separation(c2)
    PA = c1.position_angle(c2)
    if mode == "xy":
        return R.arcsec*np.sin(PA.radian), R.arcsec*np.cos(PA.radian)
    elif mode == "polar":
        return R.arcsec, PA.deg
    else:
        raise NotImplementedError

if __name__ == "__main__":
    # Coordinates of 177-341
    c1 = coord.SkyCoord("5:35:17.633", "-5:23:41.68", unit=(u.hourangle, u.deg))
    # Coordinates of th1C
    c0 = coord.SkyCoord(83.818289, -5.3895909, unit=(u.deg, u.deg))

    # Test vector separation
    x, y = vector_separation(c1, c0)
    print("Separation HST1-th1C: alpha = {} arcsec, delta = {} arcsec".format(x, y))
    # Test polar version of the same
    R, PA = vector_separation(c1, c0, mode="polar")
    print("Separation HST1-th1C: R = {} arcsec, PA = {} degrees".format(R, PA) )
    print("Scalar separation calculated by astropy: ", c1.separation(c0))

    # Coordinates of LL Ori
    c2 = coord.SkyCoord("5:35:05.566", "-5:25:19.02", unit=(u.hourangle, u.deg))
    R, PA = vector_separation(c2, c0, mode="polar")
    print("Separation LL1-th1C: R = {} arcsec, PA = {} degrees".format(R, PA) )
    print("Scalar separation calculated by astropy: ", c2.separation(c0))


test_results = """
Output generated 28 Oct 2013:

Separation proplyd-th1C: alpha = -18.5719660761 arcsec, delta = 19.15276 arcsec
Separation HST1-th1C: R = 26.6785707929 arcsec, PA = 315.8820319 degrees
Scalar separation calculated by astropy:  0d00m26.67863s
Separation LL1-th1C: R = 199.230736498 arcsec, PA = 54.2172973923 degrees
Scalar separation calculated by astropy:  0d03m19.23424s

That seems to be accurate enough
"""
