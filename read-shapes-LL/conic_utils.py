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


def yconic_th(x, th_inf=45.0):
    B = np.tan(np.radians(th_inf))
    return yconic(x, B)


def yconic(x, B=-1.0):
    """General conic section of axis ratio B=b/a"""
    if B==0.0:
        return 0.5*x**2
    else:
        s = np.sign(B)
        return s*(1.0 - np.sqrt(1.0 - s*x**2*B**2))/B**2

def x90(R0,th):
    """
    R90 value in function of R0 in the conic frame (normalized with Rc)
    Input parameters: 
    R0: Distance between the arc nose and the proplyd
    th: the assymptotic angle of hyperbolae (also applicable 
    to the other conics)
    """
    B = np.tan(np.radians(th))
    s = np.sign(B)
    if B ==0:
        x = np.sqrt(2*R0)
    else:
        x = np.sqrt((s/B**2)*( 1-B**4*( -R0+s/B**2 )**2 ) )   
    return x

    
def fit_conic(xx, yy, Rh, thh, PAh, xxh, yyh,
              symmetrical=True, freeze_theta=False, full=False, Rmin=0.0):
    """Fit conic section to the data xx, yy

    The hyperbola is described by 5 parameters:

    Rh : radius of curvature on the axis
    thh: asymptotic angle of wings from axis
    PAh: Position Angle of axis (pegged to xhh, yhh if symmetrical=True)
    xxh: xx position of center-of-curvature
    yyh: yy position of center-of-curvature

    If symmetrical=True, then the axis must pass through the star (origin of xx, yy frame)
    => tan(PAh) = xxh/yyh

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
        """Function to minimize - gives equal weight to all points and is in units of arcsec"""
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
        residual = Rh*(yconic_th(x, thh) - y)
        if full:
            return {"residual": residual, "x": x, "y": y}
        else:
            return residual

    approx_scale = abs(xx[0]) + abs(xx[-1]) + abs(yy[0]) + abs(yy[-1])

    # Pack arguments into parameters for the fit
    params = lmfit.Parameters()
    params.add("PA", value=PAh)
    params.add("R", value=Rh, min=Rmin, max=2*approx_scale)    
    params.add("xx", value=xxh)
    params.add("yy", value=yyh)
    params.add("th", value=thh, min=-90.0, max=90.0, vary=not freeze_theta)
    if symmetrical:
        params["xx"].expr = "yy * tan(PA*pi/180.0)"
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
    y = yconic_th(x, thh)
    sPA, cPA = np.sin(np.radians(PAh)), np.cos(np.radians(PAh))
    xx0, yy0 = xxh + Rh*sPA, yyh + Rh*cPA
    xx = xx0 + Rh*(x*cPA - y*sPA)
    yy = yy0 + Rh*(-x*sPA - y*cPA)
    return xx, yy

import json
import glob
import os

if __name__ == "__main__":
    print "Testing hyperbola functions ..."

    datadir = os.path.expanduser("~/Dropbox/JorgeBowshocks/")
    datafiles = glob.glob(datadir + "*/*/*-xyc.json")
    thetas = [75.0, 60.0, 45.0, 30.0, 15.0, 0.0001, -15.0, -30.0, -45.0, -60.0, -75.0, -15.0]
    lws = [1.0]*5 + [1.0] + [0.5]*5 + [2.0]
    colors = "ymcrgbgrcmyk"
    for datafile in datafiles:
        print "*"*80
        print "Processing ", datafile
        print "*"*80
        db = json.load(open(datafile))
        testdata = db["outer"]

        fig = plt.figure()
        ax = fig.add_subplot(221)
        ax.plot(testdata['x'], testdata['y'], 'ko', alpha=0.4)
        axres = fig.add_subplot(222)
        axsig = fig.add_subplot(223)
        axr90 = fig.add_subplot(224)
        model_theta_infs = []
        model_sigmas = []
        model_R90s = []
        for theta_inf, color, lw in zip(thetas, colors, lws):
            Rh, thh, PAh, xxh, yyh, fit = fit_conic(
                xx=np.array(testdata["x"]),
                yy=np.array(testdata["y"]),
                Rh=testdata["Rc"],
                thh=theta_inf,
                PAh=testdata["PAc"],
                xxh=testdata["xc"],
                yyh=testdata["yc"],
                freeze_theta=not lw==2.0,
                full=True,
                Rmin=testdata["R0"],
            )
            x, y = world_hyperbola(Rh, thh, PAh, xxh, yyh)
            R90 = Rh*x90(testdata["R0"]/Rh, thh)
            model_R90s.append(R90)
            model_theta_infs.append(thh)
            rms = np.sqrt(np.mean(fit['residual']**2))
            model_sigmas.append(rms)
            try:
                label = r"$\theta_{{\infty}} = {}^\circ$, $R_{{\mathrm{{90}}}} = {:.2f}''$".format(int(thh), R90)
                ax.plot(x, y, '-' + color, label=label, alpha=0.5, lw=lw)
                ax.arrow(xxh, yyh, 2*Rh*np.sin(np.radians(PAh)), 2*Rh*np.cos(np.radians(PAh)), 
                         fc='none', ec=color, width=0.0003, alpha=0.5, lw=0.5, head_width=0.05*Rh, head_length=0.1*Rh)
                ax.plot([xxh], [yyh], '.' + color, alpha=0.5)
                axres.plot(fit['x']*Rh, fit['residual'], ls='-', color=color, marker='.', label=label, alpha=0.5, lw=lw)
            except ValueError:
                label = "No solution"
                

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
        axres.set_xlim(None, None)
        # res_scale = 2*np.array(model_sigmas).mean()
        res_scale = testdata["R0"]
        axres.set_ylim(-res_scale, res_scale)
        axres.grid(alpha=0.2, linestyle='-')

        axsig.plot(model_theta_infs[:-1], model_sigmas[:-1], "k-")
        for th, sig, c in zip(model_theta_infs, model_sigmas, colors):
            axsig.plot(th, sig, color=c, marker="o")
        axsig.set_xlabel("theta_inf, degrees")
        axsig.set_ylabel("RMS residual, arcsec")
        axsig.set_xlim(-90.0, 90.0)
        axsig.set_ylim(0.0, None)
        axsig.grid(alpha=0.2, linestyle='-')

        axr90.plot(model_theta_infs[:-1], model_R90s[:-1], "k-")
        for th, R90, c in zip(model_theta_infs, model_R90s, colors):
            axr90.plot(th, R90, color=c, marker="o")
        axr90.plot([-90.0, 90.0], [testdata["R0"], testdata["R0"]], "k-", lw=4.0, alpha=0.3)
        axr90.set_xlabel("theta_inf, degrees")
        axr90.set_ylabel("R90, arcsec")
        axr90.set_xlim(-90.0, 90.0)
        axr90.set_ylim(0.0, None)
        axr90.grid(alpha=0.2, linestyle='-')

        fig.set_size_inches(12, 12)
        fig.tight_layout()
        figfile = os.path.basename(datafile).replace("-xyc.json", "-hyper-test.png")
        fig.savefig("test/" + figfile, dpi=300)










