import numpy as np
from scipy.integrate import odeint

def dydt(y, t, thB, LORENTZ_FAC):
    """Derivatives for ODE: x'' = 0.5 (x^{-2} - alpha^2 (x' + 1))"""
    X, U, Y, V, Z, W = y
    dXdt = U
    dYdt = V
    dZdt = W
    R2 = X**2 + Y**2 + Z**2
    # Unit vector along radius
    R_hat = np.array([X, Y, Z])/np.sqrt(R2)
    # Radiative acceleration vector
    arad = 0.5*R_hat/R2

    # B-field
    Bx, By, Bz = np.cos(thB), np.sin(thB), 0.0
    # Lorentz acceleration vector
    alorentz = 0.5*LORENTZ_FAC * np.cross([U + 1.0, V, W], [Bx, By, Bz])
    dUdt = arad[0] + alorentz[0]
    dVdt = arad[1] + alorentz[1]
    dWdt = arad[2] + alorentz[2]
    return [dXdt, dUdt, dYdt, dVdt, dZdt, dWdt]


def streamline(thB=np.radians(90), X0=10.0, Y0=0.0, Z0=0.0,
               U0=-1.0, V0=0.0, W0=0.0,
               tstop=60.0, n=201, LFAC=10.0):
    # Time grid
    t = np.linspace(0.0, tstop, n)
    # parallel stream
    # Vector of initial conditions
    y0 = [X0, U0, Y0, V0, Z0, W0]
    soln = odeint(dydt, y0, t, args=(thB, LFAC))

    return {'t': t, 'y0': Y0, 'z0': Z0, 'b': np.hypot(Y0, Z0),
            'x': soln[:, 0], 'u': soln[:, 1],
            'y': soln[:, 2], 'v': soln[:, 3],
            'z': soln[:, 4], 'w': soln[:, 5],}
