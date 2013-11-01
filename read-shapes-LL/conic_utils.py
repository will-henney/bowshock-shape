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


def fit_hyperbola(xx, yy, Rh, thh, PAh, xxh, yyh, freeze_theta=False, full=False):
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
    
    def model_minus_data(params, xx, yy, full=False):
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
        residual = yhyperbola(x, thh) - y
        if full:
            return {"residual": residual, "x": x, "y": y}
        else:
            return residual

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
    if full:
        return tuple(results + [model_minus_data(params, xx, yy, True)])
    else:
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


def axis_hyperbola(Rh, thh, PAh, xxh, yyh):
    """Anchored vector along the axis of the hyperbola

    Returns two points as [x1, x2], [y1, y2]

    Starts at the center of curvature and goes through the nose, and
    then as far again out the other side

    """
    x = [xxh, xxh + 2*Rh*np.sin(np.radians(PAh))]
    y = [yyh, yyh + 2*Rh*np.cos(np.radians(PAh))]
    return x, y


import json
import glob
import os

if __name__ == "__main__":
    print "Testing hyperbola functions ..."

    datadir = os.path.expanduser("~/Dropbox/JorgeBowshocks/")
    datafiles = glob.glob(datadir + "*/*/*-xyc.json")
    for datafile in datafiles:
        db = json.load(open(datafile))
        testdata = db["outer"]

        fig = plt.figure()
        ax = fig.add_subplot(211)
        ax.plot(testdata['x'], testdata['y'], 'ko', alpha=0.4)
        axres = fig.add_subplot(212)
        for theta_inf, color in zip([0.0001, 15.0, 30.0, 45.0, 60.0, -15.0], "bgrcmy"):
            Rh, thh, PAh, xxh, yyh, fit = fit_hyperbola(
                xx=np.array(testdata["x"]),
                yy=np.array(testdata["y"]),
                Rh=testdata["Rc"]/2,
                thh=np.abs(theta_inf),
                PAh=testdata["PAc"],
                xxh=testdata["xc"],
                yyh=testdata["yc"],
                freeze_theta=theta_inf > 0.0,
                full=True,
            )
            x, y = world_hyperbola(Rh, thh, PAh, xxh, yyh)
            label = r"$\theta_{{\infty}} = {}^\circ$, $R_{{\mathrm{{c}}}} = {:.2f}''$".format(int(thh), Rh)
            ax.plot(x, y, '-' + color, label=label, alpha=0.5)
            xxa, yya = axis_hyperbola(Rh, thh, PAh, xxh, yyh)
            ax.arrow(xxh, yyh, 2*Rh*np.sin(np.radians(PAh)), 2*Rh*np.cos(np.radians(PAh)), 
                     fc='none', ec=color, width=0.0003, alpha=0.5, lw=0.5, head_width=0.05*Rh, head_length=0.1*Rh)
            axres.plot(fit['x']*Rh, fit['residual']*Rh, '-' + color, label=label, alpha=0.5)
        ax.axis('equal')
        scale = testdata["Rc"]
        ax.set_xlim(3*scale, -3*scale)
        ax.set_ylim(-3*scale, 3*scale)
        ax.legend(loc="best", fancybox=True, shadow=True, fontsize="x-small")
        ax.set_xlabel("Delta alpha, arcsec")
        ax.set_ylabel("Delta delta, arcsec")
        ax.grid(alpha=0.2, linestyle='-')
        ax.set_title(os.path.basename(datafile))

        axres.legend(ncol=2, loc="lower center", fancybox=True, shadow=True, fontsize="x-small")
        axres.set_xlabel("Lateral coordinate, x, arcsec")
        axres.set_ylabel("(Model - data) for axial coordinate, y, arcsec")
        axres.set_xlim(-2*scale, 2*scale)
        axres.set_ylim(-0.5, 0.5)
        axres.grid(alpha=0.2, linestyle='-')

        fig.set_size_inches(6, 10)
        fig.tight_layout()
        figfile = os.path.basename(datafile).replace("-xyc.json", "-hyper-test.png")
        fig.savefig("test/" + figfile, dpi=300)










