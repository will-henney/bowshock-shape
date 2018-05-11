import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from vec_root import chandrupatla

figname = sys.argv[0].replace('.py', '.pdf')

sns.set_style('ticks')
sns.set_color_codes('bright')
fig, axes = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(8, 4))
stardata = [
    [10.0, 0.63, 0.0066, 1e-4, axes[0]],
    [20.0, 5.453, 0.1199, 0.1, axes[1]],
    [40.0, 22.19, 0.4468, 1.0, axes[2]],
]
MICRON = 1e-4                   # cm

# Velocities in units of km/s (10 km/s -> 100 km/s)
vgrid = np.linspace(10.0, 100.0, 800)
vgrid = np.logspace(1.1, 3.1, 800)
# Densities in units of 1 pcc  (0.01 -> 1e5)
logngrid = np.linspace(-3.3, 6.3, 800)
# 2d versions of velocity and density grids
vv, nn = np.meshgrid(vgrid, 10**logngrid)

def rstar(v10, n, L4):
    """Characteristic radius in pc"""
    return 2.21*np.sqrt(L4/n)/v10

def taustar(v10, n, L4, kappa600=1.0):
    """Characteristic optical depth"""
    return 0.0089*kappa600*np.sqrt(L4*n)/v10

def xfunc(x, ts, eta):
    """Function to be zeroed to find x"""
    return x**2 - (1.0 - np.exp(-2*ts*x)) - eta


def drift_velocity(G_n, Gjump=1.0e5, Gbump=10.0, wbump=2.0, wjump=10.0, jumpfac=30.0):
    """
    Fit to Cloudy model results for drift velocity as function of G/n
    in Habing cm^3.  Results in km/s 
    """
    # Linear up to the jump value
    wdrift = wjump*G_n/Gjump
    # Then much higher beyond the jump value
    wdrift[G_n > Gjump] *= jumpfac
    # Plus a bump for the low charge region
    wdrift += wbump*np.exp( -(np.log10(G_n) - np.log10(Gbump))**2 )
    return wdrift


R0s = (0.0003, 0.001, 0.003, 0.01, 0.03,
       0.1, 0.3, 1.0, 3.0, 10.0, 30.0)
lws = np.linspace(0.3, 2.5, len(R0s))
cformats = { 0.0001: "0.0001 pc", 0.001: "0.001 pc", 0.01: "0.01 pc", 0.1:
             "0.1 pc", 1.0: "1 pc", 10.0: "10 pc", }
clevs_to_label = list(cformats.keys())

box_params = dict(fc='w', ec='0.8', lw=0.4, pad=2)
RBW_label = r"Radiation bow wave, $\eta < \tau < 1$"
RBS_label = r"Radiation bow shock, $\tau > 1$"
WBS_label = r"Wind bow shock, $\tau < \eta$"

# Miscellaneous panel-dependent plot params
d = {
    "RBW y": {10.0: 4000.0, 20.0: 5000.0, 40.0: 2500.0},
    "trapped y": {10.0: 3.5e4, 20.0: 2.5e5, 40.0: 1.5e5},
    "trapped bg": {10.0: '0.85', 20.0: 'w', 40.0: 'w'},
    "IF tau": {10.0: 0.2, 20.0: 3.7, 40.0: 6.4},
    "IF tau gas": {10.0: 5.0, 20.0: 5.0, 40.0: 5.0},
}

# Data for the FUV G/n parameter
frac_fuv = {
    # Fraction of total luminosity in FUV band
    10.0: 0.718,
    20.0: 0.642,
    40.0: 0.434,
}
# Habing flux in erg/cm2/s
Habing = 1.6e-3
# Astronomical constants in cgs
Lsun = 3.82e33
pc = 3.085677582e18

# Values for graphite
Xi_dag = {
    "gra": {
        10.0: 1000.0,
        20.0: 3000.0,
        40.0: 3000.0
    },
    "sil": {
        10.0: 350.0,
        20.0: 2500.0,
        40.0: 2500.0
    },
}
Qp = 2.0
rho_d = {
    "gra": 2.2,
    "sil": 3.5,
}

T0 = 8000
kappa = 600.0
Zd = 0.01
for M, L4, eta, S49, ax in stardata:
    Mlabel = "\n".join([rf"$M = {M:.0f}\, M_\odot$",
                        rf"$L = {1e4*L4:.1e}\, L_\odot$".replace("e+0", r"\times 10^"),
                        rf"$\eta = {eta}$"])
    Rs = rstar(vv/10, nn, L4)
    ts = taustar(vv/10, nn, L4, kappa600=kappa/600.0)
    a, b = 0.0, 2*np.sqrt(1.0 + eta)
    x = chandrupatla(xfunc, a, b, args=(ts, eta))
    R0 = x*Rs
    tau = 2*x*ts

    # Radiative turnaround radius (no drag)
    Rstarstar = 2*ts*Rs / Zd

    # Ionization parameter - fiducial
    U = 2.789*S49 / (R0**2 * nn)
    # Shell ionization parameter - assume compression by M^2
    Ush = U/(vv/10)**2
    # Ionization fraction in shell
    ysh = 1.0 - 1.0/(3.5e5*Ush)


    # Fraction of ionizing photons absorbed in shell
    alphaB = 2.6e-13*(T0/1e4)**(-0.7)
    absfrac0 = 3*np.pi*(3.085677582e18)**3 * alphaB / 1e49
    absfrac = absfrac0 * (vv/10)**2 * nn**2 * R0**3 / S49
    # Equivalent optical depth
    tau_gas = -np.log(1.0 - absfrac)
    #absfrac = 2.76e-4 * (vv/10)**-1 * nn**0.5 * (L4)**1.5 / S49

    tau_dust = (3/8)*tau
    # Ionization parameter just outside the shell
    Uout = U*np.exp(-(tau_dust + tau_gas))
    Uout[~np.isfinite(Uout)] = 0.0
    Ushout = Uout/(vv/10)**2

    c_sig_over_alpha = 2.99792458e10*3e-18 / alphaB
    c_sig_over_alpha *= (1 - absfrac)**(1./3.)
    y_IF = 0.1
    y1, y2 = 0.01, 0.99
    LHS_IF = y_IF**2 / (1 - y_IF)
    LHS1, LHS2 = y1**2/(1 - y1), y2**2/(1 - y2)
    #cs = ax.contour(vv, nn, c_sig_over_alpha*Ushout, (LHS_IF,), linewidths=2, colors='r', alpha=0.5)
    m = (c_sig_over_alpha*Ushout >= LHS1) & (c_sig_over_alpha*Ushout <= LHS2)
    tau_gas_IF = np.nanmean(tau_gas[m])
    tau_dust_IF = np.nanmean(tau_dust[m])
    cs = ax.contourf(vv, nn, tau_dust, (tau_dust[m].min(), tau_dust[m].max()), linewidths=2, colors='r', alpha=0.3)
    ax.contour(vv, nn, absfrac, 1.0, linewidths=0.7, colors='r')
    ax.contourf(vv, nn, tau, (eta, 1.0), colors='k', alpha=0.2)

    for grain, a, color in [
            ["sil", 0.02, "c"],
            ["sil", 0.2, "b"],
            ["gra", 0.02, "y"],
            ["gra", 0.2, "r"],
    ]:
        kappa_d = 3.0*Qp/(4.0*a*MICRON*rho_d[grain])
        mdw = (
            (ts > 0.5*np.sqrt(eta)*kappa/kappa_d)
            & (ts < 0.5*(vv/10.0)/np.sqrt(Xi_dag[grain][M]))
            & (vv/10 > np.sqrt(eta*Xi_dag[grain][M]))
        ).astype(float)
        cs = ax.contourf(vv, nn, mdw, [0.5, 1.5], colors=color, alpha=0.4)

    cs = ax.contour(vv, nn, R0, R0s, linewidths=lws, colors='k')
    clevs = [level for level in clevs_to_label if level in cs.levels]
    ax.clabel(cs, clevs,
              fontsize='x-small', fmt=cformats,
              inline=True, inline_spacing=2, use_clabeltext=True)
    ax.text(100.0, 1e5, Mlabel, zorder=100, fontsize='x-small', bbox=box_params)

    #
    # Now do the cooling length
    #
    # Sound speed:
    cs = 11.4*np.sqrt(T0/1e4)
    # pre-shock Mach number
    M0 = vv/cs
    # post-shock Mach number


    # Now do the FUV parameter
    # cs = ax.contourf(vv, nn, G_n, (1.0e4, G_n_jump[M]),
    #                  linewidths=1, colors='g', alpha=0.15)
    # cs = ax.contourf(vv, nn, G_n, (G_n_jump[M], 1e12),
    #                 linewidths=1, colors='g', alpha=0.3)
    # cs = ax.contourf(vv, nn, G_n, (3.0, 30.0),
    #                 linewidths=1, colors='g', alpha=0.15)


    ax.set(xscale='log', yscale='log')

axes[0].set(xlabel=r"$v$, km s$^{-1}$", ylabel=r"$n$, cm$^{-3}$")
sns.despine()
for ax in axes:
    ax.label_outer()
fig.tight_layout()
fig.savefig(figname)

print(figname, end='')
