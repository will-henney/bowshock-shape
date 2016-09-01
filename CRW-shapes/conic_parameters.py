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
    """Limit of (th_1 - th_1,inf) / (th_inf - th) as th -> th_inf

    This is now also equal to J in the expansion: 

        phi_1 = J phi + K phi^2
    """
    return beta*(np.pi/G(tht) - 1.0)


def G(theta):
    """Auxiliary angle function"""
    return theta - np.sin(theta)*np.cos(theta)


def phi_ratio_anisotropic(beta, xi, tht):
    """Limit of (th_1 - th_1,inf) / (th_inf - th) as th -> th_inf"""
    # raise NotImplementedError('TODO: write phi_ratio for anisotropic case')
    return 0.5

def K_func_CRW(beta, tht, J):
    """Second order co-efficient in phi_1 = J phi + K phi^2 expansion"""
    rslt = -(1 + J)/(1 - beta)/np.tan(tht)
    rslt *= 1 - J**2*np.sin(tht)**2
    return rslt 


def K_func_anisotropic(beta, xi, tht, J):
    """Second order co-efficient in phi_1 = J phi + K phi^2 expansion"""
    # raise NotImplementedError('TODO: write K function for anisotropic case')
    return 0.0

def a_over_x(tau, J, K):
    """Hyperbola scale in terms of coefficients J and K"""
    return np.sqrt(0.5*((K/(1 + J) - tau*J)*tau/(1 + tau**2) + 4*J))


class HeadTail(object):
    """Conic fits to the head and tail"""
    def __init__(self, beta, xi=None, xmin=None):
        """`xmin` is minimum x value of natural matching point.  If `x_m` <
        `xmin`, then the matching point will be forced to be `xmin`"""
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
            # This was formerly known as phi1_over_phi
            self.J = phi_ratio_CRW(beta, self.theta_t)
            self.K = K_func_CRW(beta, self.theta_t, self.J)
            # Empirically determined correction factor
            self.K *= 0.5*self.J*(1.0 + beta)
        else:
            self.theta_t = theta_tail(beta, xi, f=finf)
            self.J = phi_ratio_anisotropic(beta, xi, self.theta_t)
            self.K = K_func_anisotropic(beta, xi, self.theta_t, self.J)
        # Center of tail hyperbola in units of R_0
        self.x0_t =  self.D / (1.0 + self.J)

        self.tau_t = np.tan(self.theta_t)
        self.T = (self.tau_h/self.tau_t)**2

        # New 30 Aug 2016
        # Scale of hyperbola now determined from the K coefficient
        # self.a_t = self.x0_t*a_over_x(self.tau_t, self.J, self.K)

        # Find the x value where two conics match in dy/dx
        self.x_m = (self.x0_t + self.sig_h*self.T*self.x0_h) / (1 + self.sig_h*self.T)
        if xmin is not None and self.x_m < xmin:
            # 30 Aug 2016: Match at x = xmin if gradients would naturally
            # match at a more negative value of x
            self.x_m = xmin
            # And throw away the previous value of x0_t so that we can
            # force y and dy/dx to match at x=0
            self.x0_t = (1 + self.sig_h*self.T)*xmin - self.sig_h*self.T*self.x0_h
      
        # Major and minor axes of tail hyperbola
        self.a_t = np.sqrt(
            (self.x_m - self.x0_t)**2
            - self.T*np.abs(self.a_h**2 - (self.x_m - self.x0_h)**2)
        )
        self.t_t = np.linspace(0.0, max(2.0, 10.0/self.a_t), 500)

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

        """Parametric Cartesian x coordinate of tail"""
        return self.x0_t - self.a_t*np.cosh(t)

    def y_tail(self, t):
        """Parametric Cartesian y coordinate of tail"""
        return self.tau_t*self.a_t*np.sinh(t)
