import sys
import json
import numpy as np
import statsmodels.api as sm
sys.path.append("../Dust-wave")
from bow_projection import Spline_R_theta_from_grid


def departure(R, theta):
    """Parabolic departure of R(theta)"""
    return 1.0/R - 0.5*(1 + np.cos(theta))

def extrapolate(mu, Delta, mu0=-0.5, force_open=False, deg=2):
    def factor(mu):
        if force_open:
            return np.abs(-1.0 - mu)**0.5
        else:
            return 1.0

    # Only fit mu < mu0
    mask = mu <= mu0
    p = np.poly1d(np.polyfit(mu[mask], Delta[mask]/factor(mu[mask]), deg=deg))
    mu_x = np.linspace(-1.0, mu0)
    return mu_x, factor(mu_x)*p(mu_x)

def R_from_Delta(mu, Delta):
    """Get radius back from departure coefficient"""
    return 1.0/(Delta + 0.5*(1.0 + mu))

class Simulation(object):
    """
    Bow shape from simulation - defined on grid and fit with splines

    Callable as function of theta
    """
    json_suffix = "-arcdata.json"
    lowess_frac = 0.2

    def __init__(self, name, npoly=2, force_open=False, mode="all"):
        self.name = name
        with open(self.name + self.json_suffix) as f:
            data = json.load(f)
        # self.thgrid = np.abs(np.radians(data['outer']['theta']))
        self.thgrid = np.radians(data['outer']['theta'])
        self.Rgrid = np.array(data['outer']['R']) / data['outer']['R0']

        if mode == "all":
            # Use all points but take absolute value of theta
            self.thgrid = np.abs(self.thgrid)
            # And do some lowess smoothing
            smooth = sm.nonparametric.lowess(self.Rgrid, self.thgrid,
                                             frac=self.lowess_frac)
            self.thgrid = smooth[:, 0]
            self.Rgrid = smooth[:, 1]
        elif mode == "positive":
            # Use only points with positive theta
            m = self.thgrid > 0.0
            self.thgrid = self.thgrid[m]
            self.Rgrid = self.Rgrid[m]
        elif mode == "negative":
            # Use only points with negative theta
            m = self.thgrid < 0.0
            self.thgrid = -self.thgrid[m]
            self.Rgrid = self.Rgrid[m]

        # Make sure arrays are sorted 
        sort_order = self.thgrid.argsort()
        self.thgrid = self.thgrid[sort_order]
        self.Rgrid = self.Rgrid[sort_order]

        # Extrapolate to theta = pi
        mu = np.cos(self.thgrid)
        Delta = departure(self.Rgrid, self.thgrid)
        mu_x, Delta_x = extrapolate(mu, Delta, deg=npoly, force_open=force_open)
        th_x = np.arccos(mu_x)
        R_x = R_from_Delta(mu_x, Delta_x)

        # Add on the extrapolated points
        self.thgrid = np.concatenate((self.thgrid, th_x))
        self.Rgrid = np.concatenate((self.Rgrid, R_x))
        # And sort again just in case
        sort_order = self.thgrid.argsort()
        self.thgrid = self.thgrid[sort_order]
        self.Rgrid = self.Rgrid[sort_order]

        # Finally do the spline fit
        self.splinefit = Spline_R_theta_from_grid(
            theta_grid=self.thgrid, R_grid=self.Rgrid)

    def __call__(self, theta):
        # When called as a function, give the spline fitted result
        return self.splinefit(theta)
