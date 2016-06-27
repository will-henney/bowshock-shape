import numpy as np
from scipy.optimize import fsolve
from scipy.special import gamma as gamma_func
from alpha_analytic import alpha_from_beta_xi as gamma

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
    numerator = (1+np.sqrt(beta))*th1_90(beta, xi)
    denominator = (1.-xi*beta)*np.sqrt(beta)
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

def finf_CRW(th, beta):
    """Function that gives f(theta) = 0 when theta = theta_infty

    Version for spherically symmetric flow, as in CRW
    """
    return th - np.tan(th) - np.pi/(1.0 -beta)


def theta_tail(beta, xi, f=finf, th_init=np.radians(91.0)):
    """Opening half-angle of tail: th1_infty

    This version only works with scalar args `beta` and `xi`
    """
    thinf = fsolve(f, th_init, args=(beta, xi))
    return np.pi - thinf
