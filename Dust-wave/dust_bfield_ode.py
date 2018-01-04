import numpy as np
from scipy.integrate import odeint

def dydt(y, t, angle):
    """Derivatives for ODE: x'' = 0.5 (x^{-2} - alpha^2 (x' + 1))"""
    X, U, Y, V = y
    dXdt = U
    dYdt = V
    R2 = X**2 + Y**2
    theta = np.arctan2(Y, X)
    dUdt = 0.0
    dVdt = 0.5*np.sin(theta)/R2
    return [dXdt, dUdt, dYdt, dVdt]


def streamline(angle=90.0, X0=10.0, Y0=0.0,
               tstop=60.0, n=201):
    # Time grid
    t = np.linspace(0.0, tstop, n)
    # parallel stream
    # Vector of initial conditions
    y0 = [X0, -1.0, Y0, 0.0]
    soln = odeint(dydt, y0, t, args=(angle,))

    return {'t': t, 'b': Y0, 
            'x': soln[:, 0], 'u': soln[:, 1],
            'y': soln[:, 2], 'v': soln[:, 3],}
