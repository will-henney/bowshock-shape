import numpy as np
from bow_projection import characteristic_radii_projected

def Rcp_R90p(inclinations, shape, *shape_args):
    """Return (~R_c', ~R_90') for `inclinations` from `shape`

`inclinations` should be a vector of angles (in radians) and `shape`
should be callable to give R(theta), optionally with additional
arguments `shape_args`

    """
    Rcp, R90p = [], []
    for inc in inclinations:
        radii = characteristic_radii_projected(inc, shape, *shape_args)
        Rcp.append(radii['tilde R_c prime'])
        R90p.append(radii['tilde R_90 prime'])
    return np.array(Rcp), np.array(R90p)


if __name__ == '__main__':
    from bow_projection import (wilkinoid_R_theta, cantoid_R_theta,
                                Spline_R_theta_from_function, theta_infinity)

    th_inf = theta_infinity(cantoid_R_theta, 0.001)
    inclinations = np.linspace(0.0, th_inf - np.pi/2, 30)
    Rc, R90 = Rcp_R90p(inclinations, cantoid_R_theta, 0.001)
    Rc_s, R90_s = Rcp_R90p(inclinations,
                           Spline_R_theta_from_function(
                               ngrid=1000,
                               shape_func=cantoid_R_theta,
                               shape_func_pars=(0.001,)))

    result = [['inc', 'R_c', 'R_c spline', 'R_90', 'R_90 spline'], None]
    result += list(zip(np.degrees(inclinations).astype(int),
                       Rc, Rc_s, R90, R90_s))
