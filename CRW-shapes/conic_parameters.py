import numpy as np
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


def B(beta, xi):
    numerator = np.sqrt(3*xi)*(1+np.sqrt(beta))
    denominator = (1.-xi*beta)*np.sqrt(1.+0.2*xi*beta)
    return numerator/denominator


def theta_c(beta,xi=1.0):
    """
    theta_c defines the excentricity of a given conic
    """
    arg = 2*A(beta,xi) - B(beta, xi)**2
    return np.sign(arg)*np.arctan(np.sqrt(np.abs(arg)))
