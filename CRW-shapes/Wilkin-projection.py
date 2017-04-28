
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import interp1d
############### Projected characteristic radii functions ###############

def f(i, Tc):
    """
    Part of equations (44), (45) and (47)
    """
    return np.sqrt(1 + Tc*np.tan(np.radians(i))**2)

def R0_proj(Rc, Tc, i):
    """
    Projected Radius at the Symmetry Axis
    """
    return 1 + Rc*(f(i, Tc) - 1)/Tc 

def elliptic_proj(Rc, Tc, i):
    """
    Projection of the ellipsoidal head. Equations (45) and (47)
    """
    Rc_p = Rc/(np.cos(np.radians(i))**2*f(i, Tc)*R0_proj(Rc, Tc, i))
    R90_p = np.sqrt(2*Rc*f(i, Tc) - Tc*R0_proj(Rc, Tc, i))/(np.cos(np.radians(i))*f(i, Tc)*np.sqrt(R0_proj(Rc, Tc, i)))
    return Rc_p, R90_p

def parabolic_proj(Rc, i):
    """
    Projection of parabolic fit to the tail. Equation (A10)
    """
    Rc_p = Rc/(np.cos(np.radians(i))**2 + 0.5*Rc*np.sin(np.radians(i))**2)
    R90_p = np.sqrt(2*Rc_p)
    return Rc_p, R90_p

######################## Real Wilkin Solution ##########################

theta = np.linspace(0, np.pi, 500, endpoint=False)
R = np.zeros_like(theta)
for j, t in enumerate(theta):
    if j == 0:
        R[j] = 1
    else:
        R[j] = np.sqrt(3*(1-t/np.tan(t))/np.sin(t)**2)

def Wilkin_R90(i, t, R):
    """
    Wilkin solution projection to obtain R'_90
    """

    def tangent_phi(x, y):
        """
        Azimutal angle tangent to LOS
        """
        tan_alpha = np.diff(y)/np.diff(x)
        return np.tan(np.radians(i))*tan_alpha
    x, y = R*np.cos(t), R*np.sin(t)
    RR, tt = R[:-1], t[:-1]
    xt = RR*(np.cos(tt)*np.cos(np.radians(i)) - np.sin(tt)*tangent_phi(x, y)*np.sin(np.radians(i)))
    yt = RR*np.sqrt(1- tangent_phi(x, y)**2)
    mask = np.isfinite(yt)
    Rt = np.hypot(xt[mask], yt[mask])
    tht = np.arctan2(yt[mask], xt[mask])
    Rt_int = interp1d(tht, Rt)
    return Rt_int(0.5*np.pi)

################ Non projected characteristic radii ####################

Rc = 5./3
R90 = np.sqrt(3.)

################ Conic parameters (Head and tail) ######################

Tch = 2*Rc - R90**2
Rct = 0.865
x0t = 2.94611

###################### Set inclination array ###########################

inc = np.linspace(0, 72)
inc_cut_1 = 55
inc_cut_2 = 70
fig = plt.figure()
ax1 = fig.add_subplot(2, 1, 1)
ax2 = fig.add_subplot(2, 1, 2)
wR90 = []
wRc = []
######################## Loop over inclination #########################
sns.set_style("white")
for i in inc:
    Rcph, R90ph = elliptic_proj(Rc, Tch, i)
    Rcpt, R90pt = parabolic_proj(Rct, i)
    ax2.plot(Rcph, Wilkin_R90(i, theta, R), "m*")
#    if i < inc_cut_1:
#        ax2.plot(Rcph, R90ph, "bo")
#    elif i < inc_cut_2:
#        ax2.plot(Rcph, R90pt, "go")
#    else:
#        ax2.plot(Rcpt, R90pt, "ro")
    wR90.append(Wilkin_R90(i, theta, R))
    wRc.append(Rcph)
Rc_grid = np.linspace(0, 10, 2000)
R90_grid = np.sqrt(2*Rc_grid)
R90_grid_s = np.sqrt(2*Rc_grid - 1.0)
R90_grid_s[~np.isfinite(R90_grid_s)] = 0.0
ax1.plot(inc, wR90, "b.")
p_fit_c = np.polyfit(np.array(np.radians(inc)), np.array(wR90), 6)
p_fit = np.poly1d(p_fit_c)
ax1.plot(np.array(inc), p_fit(np.array(np.radians(inc))), "r-", lw=2, alpha=0.5,
         label=p_fit)
ax2.plot(np.array(wRc), p_fit(np.array(np.radians(inc))))
ax2.fill_between(Rc_grid, 0, R90_grid_s, color="k", alpha=0.1)
ax2.fill_between(Rc_grid, R90_grid_s, R90_grid, color="k", alpha=0.3)
ax2.set_xlabel(r"$\tilde{R}'_c$")
ax2.set_ylabel(r"$\tilde{R}'_{90}$")
ax2.set_xlim(0, 4)
ax2.set_ylim(0, 4)
ax1.legend()
fig.savefig("Wilkin_projected.pdf")
