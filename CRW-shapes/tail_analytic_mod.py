import numpy as np

def x0_fit(beta, xi):
    """Analytic fit to x0 parameter for tail hyperbola"""
    if xi is None:
        C1 = 0.0824
        C2 = 0.4256
        C3 = 1.4298
    else:
        C1 = 0.1491 - 0.0525*xi
        C2 = 0.7779 - 0.2942*xi
        C3 = 1.9970 - 0.4385*xi

    trend = 0.7*beta**-0.55

    x = np.log10(beta)
    return trend*(C1*x**2 + C2*x + C3)


def x0_minus_a_fit(beta, xi):
    """Analytic fit to (x0 - a) parameter for tail hyperbola"""
    if xi is None:
        C1 = -0.1225
        C2 = 1.0570
    else:
        C1 = 0.0625*xi - 0.4996
        C2 = 1.1432 - 0.3691*xi
    x = np.log10(beta)
    return C1*x + C2
