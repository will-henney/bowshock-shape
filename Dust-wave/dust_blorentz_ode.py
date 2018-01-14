import numpy as np
from scipy.integrate import odeint

def dydt(y, t, thB, LORENTZ_FAC, ALPHA_DRAG=0.0):
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

    # Drag acceleration vector (assume that gas velocity is [-1, 0, 0])
    adrag = -0.5*ALPHA_DRAG**2 * np.array([U + 1.0, V, W])

    dUdt = arad[0] + alorentz[0] + adrag[0] 
    dVdt = arad[1] + alorentz[1] + adrag[1] 
    dWdt = arad[2] + alorentz[2] + adrag[2] 
    return [dXdt, dUdt, dYdt, dVdt, dZdt, dWdt]


def init_vturb(v_turb_0, mu_p_max, thB):
    """Turbulent grain velocity at a high pitch angle"""
    # Pitch angle to field
    c_p = np.random.uniform(0.0, mu_p_max)
    s_p = np.sqrt(1.0 - c_p**2) 
    # Azimuth around field
    phi_B = np.random.uniform(0.0, 2*np.pi)
    spB, cpB = np.sin(phi_B), np.cos(phi_B) 
    # Magnitude and sign
    v_turb = np.random.normal(loc=0.0, scale=v_turb_0)
    # Components in frame of B-field
    vBx = v_turb*c_p
    vBy = v_turb*s_p*cpB
    vBz = v_turb*s_p*spB
    # Now transform to flow frame
    cthB, sthB = np.cos(thB), np.sin(thB)
    u0 = vBx*cthB - vBy*sthB
    v0 = vBx*sthB + vBy*cthB
    w0 = vBz
    return u0, v0, w0


def streamline(thB=np.radians(90), X0=10.0, Y0=0.0, Z0=0.0,
               U0=-1.0, V0=0.0, W0=0.0,
               V_TURB_0=0.0, mu_p_max=0.1,
               tstop=60.0, n=201, LFAC=10.0, ALPHA_DRAG=0.0):
    # Time grid
    t = np.linspace(0.0, tstop, n)
    # Turbulent perturbation
    if V_TURB_0 > 0.0:
        u0, v0, w0 = init_vturb(V_TURB_0, mu_p_max, thB)
    else:
        u0, v0, w0 = 0.0, 0.0, 0.0  
    # Vector of initial conditions
    y0 = [X0, U0 + u0, Y0, V0 + v0, Z0, W0 + w0]
    soln = odeint(dydt, y0, t, args=(thB, LFAC, ALPHA_DRAG))

    return {'t': t, 'y0': Y0, 'z0': Z0, 'b': np.hypot(Y0, Z0),
            'x': soln[:, 0], 'u': soln[:, 1],
            'y': soln[:, 2], 'v': soln[:, 3],
            'z': soln[:, 4], 'w': soln[:, 5],}
