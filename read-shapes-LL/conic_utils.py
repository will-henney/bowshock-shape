"""
Functions for fitting conic sections to arcs

This generalises the previous circle fits

First we do hyperbolas, but maybe we will do ellipses too at a later date
"""
import numpy as np

def yhyperbol(x, th_inf=45.0):
    "Hyperbola in y(x) version"
    B = np.tan(np.radians(th_inf))
    return (np.sqrt(1.0 + x**2*B**2) - 1.0)/B**2


def yparabol(x):
    "Parabola in y(x) version"
    return 0.5*x**2


def ycircle(x):
    "Circle in y(x) version"
    return 1.0 - np.sqrt(1.0 - x**2)


def R_vs_th_hyperbola(theta, A=1.0, th_inf=45.0, B=None, units="degrees"):
    """Polar equation of a scale-free offset hyperbola 
    
    In R(theta) format with respect to an arbitrary origin
    
    A is Rc/R0, where Rc is the radius of curvature on the hyperbola
    axis and R0 is the distance of the hyperbola "nose" from the
    origin

    Returns radius in units of Rc
    """
    if B is None: 
        B = np.tan(np.radians(th_inf))

    if "deg" in units:
        th = np.radians(theta)
    else:
        th = theta

    numerator1 = (A + B**2)*np.cos(th)
    numerator2 = np.sqrt(A**2 - (A**2 + 2*A + B**2)*np.sin(th)**2)
    denominator = B**2 + (1.0 - B**2)*np.sin(th)**2
    return (numerator1 - numerator2)/denominator


if __name__ == "__main__":
    print "Testing hyperbola functions ..."
    th = np.linspace(0.0, 90.0)
    r = R_vs_th_hyperbola(th)
    for a, b in zip(th, r):
        print a, b
