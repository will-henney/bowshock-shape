import numpy as np
from scipy.optimize import brentq
from scipy.misc import derivative

# * Module parameters
#
# The delta theta that is used in the central difference approximation
# to the derivative of the R(theta) function.  For optimum balance
# between round-off and discretization error, this should be of order
# ~sqrt(eps)~, where ~eps~ is the machine precision
DX_FOR_NUMERICAL_DERIVATIVE = 3.0*np.sqrt(np.finfo(1.0).resolution)

# * Functions to find plane-of-sky shape
#
# All these functions should have argument lists of the form:
#
# :    theta, [inc], func_R, *args_for_func_R
#
# where ~func_R~ has signature ~func_R(theta, *args_for_func_R)~ and
# ~inc~ is the inclination (for those functions that depend on that).
#
# They should also be written as element-wise functions of a vector
# ~theta~, so no ~if~ statements are allowed, but ~inc~ must be a
# scalar, as must all of the extra args for ~func_R~.
#
def omega(theta, func_R, *args_for_func_R):
    """Find omega = R^{-1} d R/d theta 

Note that theta may be an array. Any extra arguments are passed to
`func_R` after `theta`

    """
    def log_R(theta, *args):
        return np.log(func_R(theta, *args))

    return derivative(log_R, theta,
                      dx=DX_FOR_NUMERICAL_DERIVATIVE, args=args_for_func_R)


def sin_phi_t(theta, inc, func_R, *args_for_func_R):
    """Returns sin(phi_t), where phi_t is azimuth along tangent line"""
    om = omega(theta, func_R, *args_for_func_R)
    tan_theta = np.tan(theta)
    return np.tan(inc)*(1.0 + om*tan_theta)/(om - tan_theta) 


def xyprime_t(theta, inc, func_R, *args_for_func_R):
    """Returns observer-frame (x', y') coordinates of tangent line"""
    R = func_R(theta, *args_for_func_R)
    sphi_t = sin_phi_t(theta, inc, func_R, *args_for_func_R)
    cos_theta, sin_theta = np.cos(theta), np.sin(theta)
    xx = cos_theta*np.cos(inc) - sin_theta*sphi_t*np.sin(inc)
    yy = sin_theta*np.sqrt(1.0 - sphi_t**2)
    return R*xx, R*yy


# * Example analytic shape functions
#

def wilkinoid_R_theta(theta):
    """Wilkin solution for wind-stream interaction"""
    return np.sqrt(3*(1.0 - theta/np.tan(theta)))/np.sin(np.abs(theta))


def cantoid_R_theta(theta, beta):
    """Cantoid solution from CRW for wind-wind interaction

Returns R(theta), normalized to the stagnation radius. Extra parameter
`beta` is the momentum ratio of the two winds

    """

    # Approximate solution for theta_1, the polar angle measured from
    # the "other" star
    theta1 = np.sqrt(7.5*(-1.0 + np.sqrt(
        1.0 + 0.8*beta*(1.0 - theta/np.tan(theta)))))

    # On-axis (theta = 0) radius to stagnation point, in units of
    # star-star separation D
    R0 = np.sqrt(beta)/(1.0+np.sqrt(beta))

    # Return radius in units of R0
    return  np.sin(theta1) / np.sin(np.abs(theta) + theta1) / R0


def paraboloid_R_theta(theta):
    """This is the special parabola with Rc = 2"""
    return 2.0 / (1.0 + np.cos(theta))


def paraboloid_omega_true(theta):
    """Analytic omega for special parabola"""
    return np.sin(theta)  / (1.0 + np.cos(theta))


# * Non-analytic shape functions
#
# These are bow shock shapes for which it is "non-trivial" to
# calculate each R(theta).  E.g., requiring numerical root finding, so
# hard to write naturally in an element-wise vector form
#
# For efficiency, we therefore calculate R(theta) once on a grid, and
# then use a spline interpolation for fast subsequent evaluation
# of R(theta) and its derivative

import scipy.interpolate

class _Spline_R_theta(object):
    """Base class for non-analytic shapes

The R(theta) shape is initialized once on a grid when the class is
instantiated, and fitted by a B-spline.  The object can then be called
as a function of theta, which will be very fast since it just
evaluates the B-spline.

    """

    thgrid = None
    Rgrid = None
    def _init_grid(self, ngrid, **shape_kwds):
        raise NotImplementedError("Override this method in a sub-class")

    def _init_spline(self, kspline, Rmax, smooth):
        """Fit a smoothing spline to the R(theta) grid. 

We fit B-splines to the parametric [x(theta), y(theta)] representation
of the bow shock. `kspline` is the order of the splines (default:
cubic). `Rmax` is the maximum radius to be included in the spline fit.
`smooth` is the spline smoothing condition (see docs for
`scipy.interpolate.splprep`).

        """
        mgood = np.isfinite(self.Rgrid) & (self.Rgrid <= Rmax)
        x = self.Rgrid[mgood]*np.cos(self.thgrid[mgood])
        y = self.Rgrid[mgood]*np.sin(self.thgrid[mgood])
        self.spline_tck, u = scipy.interpolate.splprep(
            [x, y], u=self.thgrid[mgood], k=kspline, s=smooth)

    def __call__(self, theta):
        """Evaluate R(theta) from spline interpolation"""
        x, y = scipy.interpolate.splev(theta, self.spline_tck)
        return np.hypot(x, y)

    def __init__(self, ngrid=100, kspline=3, Rmax=100, smooth=0, **shape_kwds):
        """"""
        # Set up grids of theta and R
        self._init_grid(ngrid, **shape_kwds)
        # Set up spline interpolation
        self._init_spline(kspline, Rmax, smooth)


class Spline_R_theta_from_function(_Spline_R_theta):
    """Spline-interpolated bow shock shape from explicit function

Extra parameters for initialization: `shape_func` and
`shape_func_pars`. THIS IS FOR TESTING ONLY!!! It checks that the
interpolation machinery works for simple shapes. Outside of such
tests, there is really no need to use the spline interpolation for
these cases.

    """

    def _init_grid(self, ngrid,
                   shape_func=paraboloid_R_theta,
                   shape_func_pars=()):
        # Include the negative branch so the spline will have the
        # right gradient on the axis
        self.thgrid = np.linspace(-np.pi, np.pi, ngrid)
        self.Rgrid = shape_func(self.thgrid, *shape_func_pars)


class Spline_R_theta_from_grid(_Spline_R_theta):
    """Spline-interpolated bow shock shape from user-specified arrays

Extra parameters for initialization: `theta_grid` and `R_grid`

    """
    def _init_grid(self, ngrid, theta_grid=None, R_grid=None):
        # Note that ngrid is ignored in this implementation
        if theta_grid is not None and R_grid is not None:
            self.thgrid = theta_grid
            self.Rgrid = R_grid
        else:
            raise ValueError("Both theta_grid and R_grid must be specified")


# * Basic tests of functionality
#

if __name__ == "__main__":
    import sys
    from matplotlib import pyplot as plt
    import seaborn as sns

    lib_name = sys.argv[0].replace('.py', '')

    sns.set_style('ticks')
    fig, ax = plt.subplots()

    th = np.linspace(0.0, np.pi, 501)
    th_dg = np.degrees(th)
    ax.plot(th_dg, omega(th, paraboloid_R_theta),
            label="paraboloid")
    ax.plot(th_dg, omega(th, wilkinoid_R_theta),
            label="wilkinoid")
    for beta in 0.001, 0.01, 0.1:
        ax.plot(th_dg, omega(th, cantoid_R_theta, beta),
                label=fr"cantoid $\beta = {beta:.3f}$")
    ax.legend(title=r"Analytic $R(\theta)$ functions")
    ax.axhline(1.0, xmin=0.35, xmax=0.65, color='white', lw=4, zorder=100)
    ax.axhline(1.0, xmin=0.35, xmax=0.65, color='k', lw=1, ls=':', zorder=101)
    ax.axhspan(0.0, 1.0, color='k', alpha=0.05, ec='none')
    ax.set_yscale('symlog', linthreshy=1.0, linscaley=0.5)
    ax.set(
        xlabel=r"Polar angle: $\theta$, degrees",
        ylabel=r"$\omega \equiv R^{-1} d R / d \theta$",
        xlim=[0, 180],
        ylim=[0.0, 2e2],
        xticks=[0, 30, 60, 90, 120, 150, 180],
    )
    sns.despine()
    fig.tight_layout()
    figfile = f"test_{lib_name}_omega.pdf"
    fig.savefig(figfile)
    print(figfile, end='')
