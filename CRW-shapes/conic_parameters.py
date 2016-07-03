import sys
import numpy as np
from scipy.optimize import fsolve
from scipy.special import gamma as gamma_func
from alpha_analytic import alpha_from_beta_xi as gamma
sys.path.append('../conic-projection')
from conproj_utils import Conic

def A(beta, xi):
    return 1.0/(1.0 - 2.0*gamma(beta, xi))


def th1_90_method1(beta, xi):
    x = 3*beta*xi
    return np.sqrt(x/(1 + x/5))


def th1_90_method2(beta, xi):
    return np.sqrt(2.5*(np.sqrt(1.0 + 12*xi*beta/5.0) - 1.0))


def B(beta,xi=1.0, th1_90=th1_90_method1):
    """
    Returns R_90 normalized with R_0
    """
    if xi is None:
        xxi = 1.0
    else:
        xxi = xi
    numerator = (1+np.sqrt(beta))*th1_90(beta, xxi)
    denominator = (1.-xxi*beta)*np.sqrt(beta)
    return numerator/denominator


def theta_c(beta,xi=1.0):
    """
    theta_c defines the excentricity of a given conic
    """
    arg = 2*A(beta,xi) - B(beta, xi)**2
    return np.sign(arg)*np.arctan(np.sqrt(np.abs(arg)))



#
# Now, functions for the hyperbola that fits the tail
#

def finf(th, beta, xi):
    """Function that gives f(theta) = 0 when theta = theta_infty

    Version for hemispheric flow with anisotropy xi
    """
    k = 2./xi-2
    C = (k+2*(1-beta))/(k+2)
    I = np.sqrt(np.pi)*gamma_func(0.5*(k+1))/(4*gamma_func(0.5*k+2))
    D = np.pi + 2*beta*I
    return th - C*np.tan(th) - D

def finf_CRW(th, beta, xi):
    """Function that gives f(theta) = 0 when theta = theta_infty

    Version for spherically symmetric flow, as in CRW
    """
    assert xi is None, 'Parameter xi is meaningless for vanilla CRW'
    return th - np.tan(th) - np.pi/(1.0 -beta)


def theta_tail(beta, xi, f=finf, th_init=np.radians(91.0)):
    """Opening half-angle of tail: th1_infty

    This version only works with scalar args `beta` and `xi`
    """
    thinf, = fsolve(f, th_init, args=(beta, xi))
    return np.pi - thinf


def phi_ratio_CRW(beta, tht):
    """Limit of (th_1 - th_1,inf) / (th_inf - th) as th -> th_inf"""
    return beta*(np.pi/(tht - np.sin(tht)*np.cos(tht)) - 1.0)


def phi_ratio_anisotropic(beta, xi, tht):
    """Limit of (th_1 - th_1,inf) / (th_inf - th) as th -> th_inf"""
    raise NotImplementedError('TODO: write phi_ratio for anisotropic case')


class HeadTail(object):
    """Conic fits to the head and tail"""
    def __init__(self, beta, xi=None):
        #
        # Set head parameters
        #
        self.A_h = A(beta, xi)
        self.theta_h = theta_c(beta, xi)
        self.sig_h = np.sign(self.theta_h)
        self.tau_h = np.tan(np.abs(self.theta_h))
        self.a_h = self.A_h/self.tau_h**2
        self.x0_h = 1.0 - self.sig_h*self.a_h
        self.head_conic = Conic(A=self.A_h, th_conic=np.degrees(self.theta_h))
        self.t_h = self.head_conic.make_t_array()
        #
        # Set tail parameters
        #
        # Distance to other star in units of R0
        self.D = (1 + np.sqrt(beta)) / np.sqrt(beta)
        if xi == None:
            # Opening angle of tail
            self.theta_t = theta_tail(beta, xi, f=finf_CRW)
            # Center of tail hyperbola in units of R_0
            self.phi1_over_phi = phi_ratio_CRW(beta, self.theta_t)
        else:
            self.theta_t = theta_tail(beta, xi, f=finf)
            self.phi1_over_phi = phi_ratio_anisotropic(beta, xi, self.theta_t)
        self.x0_t =  self.D / (1.0 + self.phi1_over_phi)
        self.tau_t = np.tan(self.theta_t)
        self.T = (self.tau_h/self.tau_t)**2
        # Find the x value where two conics match in y and dy/dx
        self.x_m = (self.x0_t + self.sig_h*self.T*self.x0_h) / (1 + self.sig_h*self.T)
        # Major and minor axes of tail hyperbola
        self.a_t = np.sqrt(self.T*self.a_h**2
                           + self.T*(self.T - 1)*(self.x_m - self.x0_h)**2)
        self.t_t = np.linspace(0.0, 10.0, 100)

    def x_head(self, t):
        """Parametric Cartesian x coordinate of head"""
        return self.head_conic.x(t)

    def y_head(self, t):
        """Parametric Cartesian y coordinate of tail"""
        return self.head_conic.y(t)

    def x_tail(self, t):
        """Parametric Cartesian x coordinate of tail"""
        return self.x0_t - self.a_t*np.cosh(t)

    def y_tail(self, t):
        """Parametric Cartesian y coordinate of tail"""
        return self.tau_t*self.a_t*np.sinh(t)
