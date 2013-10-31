"""Functions for fitting conic sections to arcs

This generalises the previous circle fits

First we do hyperbolas, but maybe we will do ellipses too at a later date

CONVENTIONS USED IN THIS FILE 
-----------------------------

* (xx, yy) are Cartesian world coordinates with xx increasing to left
  (for instance RA and Dec arcsec offsets from star)

* (x, y) are normalized Cartesian coordinates in the frame of the
  hyperbola (or other conic section), with x increasing to right

Note that different conventions are used elsewhere

"""

import numpy as np
import lmfit
import matplotlib.pyplot as plt


def yhyperbola(x, th_inf=45.0):
    "Hyperbola in y(x) version"
    B = np.tan(np.radians(th_inf))
    return (np.sqrt(1.0 + x**2*B**2) - 1.0)/B**2


def yparabola(x):
    "Parabola in y(x) version"
    return 0.5*x**2


def ycircle(x):
    "Circle in y(x) version"
    return 1.0 - np.sqrt(1.0 - x**2)


def fit_hyperbola(xx, yy, Rh, thh, PAh, xxh, yyh, freeze_theta=False):
    """Fit hyperbola to the data xx, yy

    The hyperbola is described by 5 parameters:

    Rh : radius of curvature on the axis
    thh: asymptotic angle of wings from axis
    PAh: Position Angle of axis
    xxh: xx position of center-of-curvature
    yyh: yy position of center-of-curvature

    The hyperbola is given as y(x), where y goes along the hyperbola
    axis and x is perpendicular to the hyperbola axis.

    The origin of the xy frame (x = y = 0) is the "nose" of the
    hyperbola, with coordinates in the xx-yy frame of 

    xx0 = xxh + Rh sin(PAh)
    yy0 = yyh + Rh cos(PAh)

    The y-axis is along the [xx, yy] unit vector: 

    yhat = [-sin(PAh), -cos(PAh)]

    The x-axis is along the [xx, yy] unit vector: 

    xhat = [-cos(PAh), -sin(PAh)]

    Unit distance in the xy frame corresponds to Rh in the xx-yy frame

    So that x = ((xx - xx0) (-cos(PAh)) + (yy - yy0) (-sin(PAh)))/Rh
            y = ((xx - xx0) (-sin(PAh)) + (yy - yy0) (-cos(PAh)))/Rh

    """
    
    def model_minus_data(params, xx, yy):
        """Function to minimize - gives equal weight to all points"""
        # Unpack parameters
        PAh = params["PA"].value
        Rh = params["R"].value
        xxh = params["xx"].value
        yyh = params["yy"].value
        thh = params["th"].value
        # Transform from (xx,yy) to (x,y) frame
        sPA, cPA = np.sin(np.radians(PAh)), np.cos(np.radians(PAh))
        xx0, yy0 = xxh + Rh*sPA, yyh + Rh*cPA
        x = (-(xx0 - xx)*cPA + (yy0 - yy)*sPA)/Rh
        y = ((xx0 - xx)*sPA + (yy0 - yy)*cPA)/Rh
        return yhyperbola(x, thh) - y

    # Pack arguments into parameters for the fit
    params = lmfit.Parameters()
    params.add("PA", value=PAh)
    params.add("R", value=Rh)
    params.add("xx", value=xxh)
    params.add("yy", value=yyh)
    params.add("th", value=thh, min=0.0, max=90.0, vary=not freeze_theta)
    lmfit.minimize(model_minus_data, params, args=(xx, yy))
    lmfit.report_errors(params)

    # Unpack parameters again for results to return
    results = [params[k].value for k in "R", "th", "PA", "xx", "yy"]
    return tuple(results)


def world_hyperbola(Rh, thh, PAh, xxh, yyh, xmin=-2.0, xmax=2.0, N=200):
    """Return the (xx, yy) world coordinates of a hyperbola

    Required arguments as in fit_hyperbola: 
    Rh - radius of curvature
    thh - asymptotic angle
    PAh - orientation of axis
    (xxh, yyh) - location of center of curvature
    """
    x = np.linspace(xmin, xmax, N) 
    y = yhyperbola(x, thh)
    sPA, cPA = np.sin(np.radians(PAh)), np.cos(np.radians(PAh))
    xx0, yy0 = xxh + Rh*sPA, yyh + Rh*cPA
    xx = xx0 + Rh*(x*cPA - y*sPA)
    yy = yy0 + Rh*(-x*sPA - y*cPA)
    return xx, yy


testdata = {
    "xc": -5.665723355313418, 
    "yc": -0.8801883891895861, 
    "Rc": 8.720647445668673, 
    "PAc": 81.16950464679529, 
    "x": [
        -1.6874162965569341, -0.9258390297935146, -0.2687919763362173, 
        0.4629195150877455, 1.0453021307693653, 1.5679531959630597, 
        2.060738486199274, 2.6431211016899057, 2.9716446285140483, 
        3.0761748416291823, 3.0313761788109854, 2.867114415398914, 
        2.672986876838374, 2.404194900502157, 1.9562082730841397, 
        1.2394296693299052, 0.6122483908300884, -0.16426176322108316, 
        -1.0602350180571172
    ], 
    "y": [
        7.310000000000727, 6.820000000001514, 6.050000000000466, 
        5.310000000000414, 4.620000000003088, 3.830000000000311, 
        3.090000000000259, 1.9900000000010465, 0.8300000000030394, 
        -0.2299999999991087, -1.2799999999987932, -2.2799999999989495, 
        -3.2299999999995777, -4.2999999999977945, -5.179999999997165, 
        -6.179999999997321, -7.029999999998893, -7.959999999997791, 
        -8.749999999997371
    ], 
}


if __name__ == "__main__":
    print "Testing hyperbola functions ..."
    plt.plot(testdata['x'], testdata['y'], 'o')

    for theta_inf in 0.0001, 15.0, 30.0, 45.0:
        Rh, thh, PAh, xxh, yyh = fit_hyperbola(
            xx=np.array(testdata["x"]),
            yy=np.array(testdata["y"]),
            Rh=testdata["Rc"]/3,
            thh=theta_inf,
            PAh=testdata["PAc"],
            xxh=testdata["xc"],
            yyh=testdata["yc"],
            freeze_theta=True
        )
        x, y = world_hyperbola(Rh, thh, PAh, xxh, yyh)
        plt.plot(x, y, '-', label=str(int(thh)))
    plt.axis('equal')
    plt.xlim(10, -10)
    plt.ylim(-10, 10)
    plt.legend()
    plt.savefig("test/hyperbola-test.pdf")

    









