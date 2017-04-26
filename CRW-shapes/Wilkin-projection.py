
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

############### Projected characteristic radii functions ###############

def f(i, Tc):
    """
    Part of equations (44), (45) and (47)
    """
    return np.sqrt(1 + Tc*np.tan(np.radians(i)))

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
######################## Loop over inclination #########################
sns.set_style("white")
for i in inc:
    Rcph, R90ph = elliptic_proj(Rc, Tch, i)
    Rcpt, R90pt = parabolic_proj(Rct, i)

    if i < inc_cut_1:
        plt.plot(Rcph, R90ph, "bo")
    elif i < inc_cut_2:
        plt.plot(Rcph, R90pt, "go")
    else:
        plt.plot(Rcpt, R90pt, "ro")

Rc_grid = np.linspace(0, 10, 2000)
R90_grid = np.sqrt(2*Rc_grid)
R90_grid_s = np.sqrt(2*Rc_grid - 1.0)
R90_grid_s[~np.isfinite(R90_grid_s)] = 0.0
plt.fill_between(Rc_grid, 0, R90_grid_s, color="k", alpha=0.1)
plt.fill_between(Rc_grid, R90_grid_s, R90_grid, color="k", alpha=0.3)
plt.xlabel(r"$\tilde{R}'_c$")
plt.ylabel(r"$\tilde{R}'_{90}$")
plt.xlim(0, 4)
plt.ylim(0, 4)
plt.savefig("Wilkin_projected.pdf")
