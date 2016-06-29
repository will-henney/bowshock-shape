"""Miscellaneous functions for Canto, Raga, Wilkin (1996) bow shocks

"""

import numpy as np


def th1_approx(th, beta):
    """Equation (26) of CRW"""
    fac = 0.8*beta*(1.0 - th/np.tan(th))
    return np.sqrt(7.5*(np.sqrt(1.0 + fac) - 1.0))
  

def radius(th, th1):
    """Radius in terms of D from Eq (23) of CRW
  
    This applies generally, even for anisotropic cases
    """
    return np.sin(th1)/np.sin(th + th1)
