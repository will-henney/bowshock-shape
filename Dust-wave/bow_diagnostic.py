import numpy as np
from astropy.table import Table
from bow_projection import characteristic_radii_projected

def parameter_table(inclinations, shape, *shape_args):
    """Diagnostic parameters for `inclinations` from `shape`

Input argument `inclinations` should be a vector of angles (in
radians) and input argument `shape` should be callable to give
R(theta), optionally with additional arguments `shape_args`.  Example
functions suitable for passing as `shape` can be found in the
`bow_projection` module.

Returns an `astropy.table.Table` of characteristic angles and radii.

    """
    rows = [characteristic_radii_projected(inc, shape, *shape_args)
            for inc in inclinations]
    tab = Table(rows=rows)
    tab['inc'] = inclinations
    return tab


if __name__ == '__main__':
    from bow_projection import (wilkinoid_R_theta, cantoid_R_theta,
                                Spline_R_theta_from_function, theta_infinity)

    th_inf = theta_infinity(cantoid_R_theta, 0.001)
    inclinations = np.linspace(0.0, th_inf - np.pi/2, 30)
    tab = parameter_table(inclinations, cantoid_R_theta, 0.001)
    Rc, R90 = tab['tilde R_c prime'], tab['tilde R_90 prime']
    tab_s = parameter_table(inclinations,
                            Spline_R_theta_from_function(
                                ngrid=1000,
                                shape_func=cantoid_R_theta,
                                shape_func_pars=(0.001,)))
    Rc_s, R90_s = tab_s['tilde R_c prime'], tab_s['tilde R_90 prime']

    result = [['inc', 'R_c', 'R_c spline', 'R_90', 'R_90 spline'], None]
    result += list(zip(np.degrees(inclinations).astype(int),
                       Rc, Rc_s, R90, R90_s))
