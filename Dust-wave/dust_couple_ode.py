import numpy as np
from scipy.integrate import odeint

def dydt(y, t, alpha):
    """Derivatives for ODE: x'' = 0.5 (x^{-2} - alpha^2 (x' + 1))"""
    X, U, Y, V = y
    dXdt = U
    dYdt = V
    R2 = X**2 + Y**2
    theta = np.arctan2(Y, X)
    dUdt = 0.5*(np.cos(theta)/R2 - alpha**2 * (U + 1.0))
    dVdt = 0.5*(np.sin(theta)/R2 - alpha**2 * V)
    return [dXdt, dUdt, dYdt, dVdt]


def streamline(alpha=1.0/3.0, x0=10.0, y0=0.0,
               tstop=60.0, n=201):
    # Initial conditions
    y0 = [x0, -1.0, y0, 0.0]
    # Time grid
    t = np.linspace(0.0, tstop, n)
    soln = odeint(dydt, y0, t, args=(alpha,))
    return {'t': t, 'b': y0, 'alpha': alpha,
            'x': soln[:, 0], 'u': soln[:, 1],
            'y': soln[:, 2], 'v': soln[:, 3],}
