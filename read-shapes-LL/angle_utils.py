"""
Functions to make it easier to use astropy.coordinates objects
"""
from astropy import coordinates as coord
from astropy import units as u
import numpy as np

def vector_separation(c1, c2, mode="xy"):
    """
    Find vector separation between two coordinates c1, c2

    Input arguments c1, c2 should be astropy.coordinates objects (or similar).
    When mode="xy" (default), return the RA and Dec separations in arcsec.
    When mode="polar", return the scalar separation in arcsec and PA of the separation.
    In either case, pure numbers are returned.  
    The sense of the separation is that it takes you from c1 to c2.
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

if __name__ == "__main__":
    # Coordinates of 177-341
    c1 = coord.ICRSCoordinates("5:35:17.633", "-5:23:41.68", unit=(u.hour, u.deg))
    # Coordinates of th1C
    c2 = coord.ICRSCoordinates(83.818289, -5.3895909, unit=(u.deg, u.deg))

    # Test vector separation
    x, y = vector_separation(c1, c2)
    print "Separation proplyd-th1C: alpha = {} arcsec, delta = {} arcsec".format(x, y) 
    # Test polar version of the same
    R, PA = vector_separation(c1, c2, mode="polar")
    print "Separation proplyd-th1C: R = {} arcsec, PA = {} degrees".format(R, PA) 

    print "Scalar separation calculated by astropy: ", c1.separation(c2)


test_resuts = """
Output generated 27 Oct 2013:
Separation proplyd-th1C: alpha = -18.5719660761 arcsec, delta = 19.15276 arcsec
Separation proplyd-th1C: R = 26.6785707929 arcsec, PA = 315.8820319 degrees
Scalar separation calculated by astropy:  0d00m26.67863s

That seems to be accurate enough
"""
