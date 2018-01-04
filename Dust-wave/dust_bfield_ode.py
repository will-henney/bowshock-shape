import numpy as np
from scipy.integrate import odeint

def dydt(y, t, thB):
    """Derivatives for ODE: x'' = 0.5 (x^{-2} - alpha^2 (x' + 1))"""
    X, U, Y, V, Z, W = y
    dXdt = U
    dYdt = V
    dZdt = W
    R2 = X**2 + Y**2 + Z**2
    # Radiative acceleration
    arad = 0.5/R2
    # dot product of radial vector with b vector
    cos_delta = (X*np.cos(thB) + Y*np.sin(thB))/np.sqrt(R2)
    # Projection onto B direction
    a_para = arad*cos_delta
    dUdt = a_para*np.cos(thB)
    dVdt = a_para*np.sin(thB)
    dWdt = 0.0
    return [dXdt, dUdt, dYdt, dVdt, dZdt, dWdt]


def streamline(thB=np.radians(90), X0=10.0, Y0=0.0, Z0=0.0,
               tstop=60.0, n=201):
    # Time grid
    t = np.linspace(0.0, tstop, n)
    # parallel stream
    # Vector of initial conditions
    y0 = [X0, -1.0, Y0, 0.0, Z0, 0.0]
    soln = odeint(dydt, y0, t, args=(thB,))

    return {'t': t, 'y0': Y0, 'z0': Z0, 'b': np.hypot(Y0, Z0),
            'x': soln[:, 0], 'u': soln[:, 1],
            'y': soln[:, 2], 'v': soln[:, 3],
            'z': soln[:, 4], 'w': soln[:, 5],}
