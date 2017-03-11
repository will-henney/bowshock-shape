import numpy as np
import astropy.coordinates as coord
import astropy.units as u
import lmfit

def Rc_from_data(points, center):
    return np.mean(center.separation(points))

def deviation_from_circle(points, center):
    return center.separation(points) - Rc_from_data(points, center)

def model_minus_data(params, points):
    center = coord.SkyCoord(params["ra"].value*u.deg, params["dec"].value*u.deg)
    return deviation_from_circle(points, center).arcsec

def fit_circle(points, center0):
    """Fit a circle to `points` with initial guess that center is at
`center0`.  Returns radius of curvature and center of curvature"""
    params = lmfit.Parameters()
    params.add("ra", value=center0.ra.deg)
    params.add("dec", value=center0.dec.deg)
    out = lmfit.minimize(model_minus_data, params, args=(points,))
    lmfit.report_fit(out)
    center = coord.SkyCoord(out.params["ra"].value*u.deg, out.params["dec"].value*u.deg)
    Rc = Rc_from_data(points, center)
    return Rc, center
