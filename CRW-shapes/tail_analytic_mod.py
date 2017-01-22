import numpy as np
from astropy.table import Table

xitab = Table.read('params-tab.tab', format='ascii.tab')
isotab = Table.read('params-iso-tab.tab', format='ascii.tab')

def x0_fit(beta, xi):
    """Analytic fit to x0 parameter for tail hyperbola"""
    if xi is None:
        pb = np.poly1d(isotab['x0_coef'])
    else:
        C3 = np.poly1d(xitab['3_ord_x0'])(xi)
        C2 = np.poly1d(xitab['2_ord_x0'])(xi)
        C1 = np.poly1d(xitab['1_ord_x0'])(xi)
        C0 = np.poly1d(xitab['0_ord_x0'])(xi)
        pb = np.poly1d([C3, C2, C1, C0])
    trend = 0.7*beta**-0.55
    x = np.log10(beta)
    return trend*pb(x)


def x0_minus_a_fit(beta, xi):
    """Analytic fit to (x0 - a) parameter for tail hyperbola"""
    if xi is None:
        pb = np.poly1d(isotab['x0Ma_coef'])
    else:
        C2 = np.poly1d(xitab['2_ord_x0Ma'])(xi)
        C1 = np.poly1d(xitab['1_ord_xoMa'])(xi)
        C0 = np.poly1d(xitab['0_ord_x0Ma'])(xi)
        pb = np.poly1d([C2, C1, C0])
    x = np.log10(beta)
    return pb(x)
