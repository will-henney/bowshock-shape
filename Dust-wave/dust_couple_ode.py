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

def dydt_div(y, t, alpha, mu):
    """Derivatives for ODE in divergent case"""
    X, U, Y, V = y
    dXdt = U
    dYdt = V
    R2 = X**2 + Y**2
    # Gas flow diverges from point (X1, Y1) = (1/mu, 0)
    X1, Y1 = 1.0/mu, 0.0
    # Gas flow is radial from that point
    R1 = np.hypot(X - X1, Y - Y1)
    U1 = (X - X1)/R1
    V1 = (Y - Y1)/R1
    theta = np.arctan2(Y, X)
    dUdt = 0.5*(np.cos(theta)/R2 - alpha**2 * (U - U1))
    dVdt = 0.5*(np.sin(theta)/R2 - alpha**2 * (V - V1))
    return [dXdt, dUdt, dYdt, dVdt]


def streamline(alpha=1.0/3.0, X0=10.0, Y0=0.0,
               tstop=60.0, n=201, mu=None):
    # Time grid
    t = np.linspace(0.0, tstop, n)
    if mu is None:
        # parallel stream
        # Vector of initial conditions
        y0 = [X0, -1.0, Y0, 0.0]
        soln = odeint(dydt, y0, t, args=(alpha,))
    else:
        # divergent stream
        X1, Y1 = 1.0/mu, 0.0
        assert X0 < X1, 'Start point must be to left of wind source'
        R1 = np.hypot(X0 - X1, Y0 - Y1)
        U0 = (X0 - X1)/R1
        V0 = (Y0 - Y1)/R1
        # Vector of initial conditions
        y0 = [X0, U0, Y0, V0]
        soln = odeint(dydt_div, y0, t, args=(alpha, mu))

    return {'t': t, 'b': Y0, 'alpha': alpha, 'mu': mu, 
            'x': soln[:, 0], 'u': soln[:, 1],
            'y': soln[:, 2], 'v': soln[:, 3],}
