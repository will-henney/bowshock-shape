import sys
import json
import numpy as np
import statsmodels.api as sm
from astropy.coordinates import Longitude
from astropy.modeling import models, fitting
sys.path.append("../Dust-wave")
from bow_projection import Spline_R_theta_from_grid


def departure(R, theta):
    """Parabolic departure of R(theta)"""
    return 1.0/R - 0.5*(1 + np.cos(theta))


def R_from_Delta(mu, Delta):
    """Get radius back from departure coefficient"""
    return 1.0/(Delta + 0.5*(1.0 + mu))


JSON_SUFFIX = "-arcdata.json"

def load_R_th(arc_prefix):
    jfile = arc_prefix + JSON_SUFFIX
    data = json.load(open(jfile))
    R0 = np.array(data['outer']['R0'])
    R = np.array(data['outer']['R'])
    th = np.radians(data['outer']['theta'])
#    th = Longitude(data['outer']['theta'], unit='deg')
#    th += Longitude(data['outer']['PA0'], unit='deg')
    return th, R/R0


class Simulation(object):
    """
    Bow shape from simulation - defined on grid and fit with splines

    Callable as function of theta
    """
    lowess_frac = 0.2

    def extrapolation_factor(self, mu):
        if self.force_open:
            return np.abs(-1.0 - mu)**0.5
        else:
            return 1.0

    def extrapolation(self, mu):
        return self.extrapolation_factor(mu)*self.extrap_polyfit(mu)

    def __init__(self, name, extrap_degree=2, mu0=-0.5,
                 cheby_degree=10, force_open=False, mode="all"):
        self.name = name
        self.force_open = force_open
        self.thgrid, self.Rgrid = load_R_th(name)
        self.thmax = self.thgrid.max()

        # Set up grid of departure function vs mu
        Delta = departure(self.Rgrid, self.thgrid)
        mu = np.cos(self.thgrid)

        # Set up Chebyshev fit to grid data (theta < thmax)
        self.chebyfit = models.Chebyshev1D(degree=cheby_degree)
        fitter = fitting.LevMarLSQFitter()
        self.chebyfit = fitter(self.chebyfit, mu, Delta)

        # Set up extrapolation fit for theta > thmax
        # Only fit mu < mu0
        mask = mu <= mu0
        self.extrap_polyfit = np.poly1d(np.polyfit(
            mu[mask], Delta[mask]/self.extrapolation_factor(mu[mask]),
            deg=extrap_degree))

        # if mode == "all":
        #     # Use all points but take absolute value of theta
        #     self.thgrid = np.abs(self.thgrid)
        #     # And do some lowess smoothing
        #     smooth = sm.nonparametric.lowess(self.Rgrid, self.thgrid,
        #                                      frac=self.lowess_frac)
        #     self.thgrid = smooth[:, 0]
        #     self.Rgrid = smooth[:, 1]
        # elif mode == "positive":
        #     # Use only points with positive theta
        #     m = self.thgrid > 0.0
        #     self.thgrid = self.thgrid[m]
        #     self.Rgrid = self.Rgrid[m]
        # elif mode == "negative":
        #     # Use only points with negative theta
        #     m = self.thgrid < 0.0
        #     self.thgrid = -self.thgrid[m]
        #     self.Rgrid = self.Rgrid[m]

        # # Make sure arrays are sorted 
        # sort_order = self.thgrid.argsort()
        # self.thgrid = self.thgrid[sort_order]
        # self.Rgrid = self.Rgrid[sort_order]

        # th_x = np.arccos(mu_x)
        # R_x = R_from_Delta(mu_x, Delta_x)

        # # Add on the extrapolated points
        # self.thgrid = np.concatenate((self.thgrid, th_x))
        # self.Rgrid = np.concatenate((self.Rgrid, R_x))
        # # And sort again just in case
        # sort_order = self.thgrid.argsort()
        # self.thgrid = self.thgrid[sort_order]
        # self.Rgrid = self.Rgrid[sort_order]

        # # Finally do the spline fit
        # self.splinefit = Spline_R_theta_from_grid(
        #     theta_grid=self.thgrid, R_grid=self.Rgrid)

    def __call__(self, theta):
        # When called as a function, give the fitted result
        mu = np.cos(theta)
        # Use Chebyshev for the range of the grid data
        # and use extrapolation for larger angles
        mask = np.cos(theta) >= np.cos(self.thmax)
        Delta = np.empty_like(mu)
        Delta[mask] = self.chebyfit(mu[mask])
        Delta[~mask] = self.extrapolation(mu[~mask])
        return R_from_Delta(mu, Delta)
