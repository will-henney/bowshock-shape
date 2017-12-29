import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from vec_root import chandrupatla

figname = sys.argv[0].replace('.py', '.pdf')

sns.set_style('ticks')
sns.set_color_codes('dark')
fig, axes = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(3.3, 4))
stardata = [
    [20.0, 15.6, 0.0476, 0.0, axes],
]

# Velocities in units of km/s (10 km/s -> 100 km/s)
vgrid = np.linspace(10.0, 100.0, 800)
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
    "RBW y": {33.0: 1500.0, 20.0: 1e5, 40.0: 2500.0},
    "trapped y": {33.0: 250, 20.0: 2.5e5, 40.0: 1.5e5},
    "trapped bg": {33.0: 'w', 20.0: 'w', 40.0: 'w'},
    "IF tau": {33.0: 0.2, 20.0: 3.7, 40.0: 6.4},
    "IF tau gas": {33.0: 5.0, 20.0: 5.0, 40.0: 5.0},
}

T0 = 1000
kappa = 60.0
for M, L4, eta, S49, ax in stardata:
    Mlabel = "\n".join([
        "M supergiant", "",
        rf"$M = {M:.0f}\, M_\odot$",
        rf"$L = {1e4*L4:.1e}\, L_\odot$".replace("e+0", r"\times 10^"),
        rf"$\eta = {eta}$"])
    Rs = rstar(vv/10, nn, L4)
    ts = taustar(vv/10, nn, L4, kappa600=kappa/600.0)
    a, b = 0.0, 2*np.sqrt(1.0 + eta)
    x = chandrupatla(xfunc, a, b, args=(ts, eta))
    R0 = x*Rs
    tau = 2*x*ts


    ax.contourf(vv, nn, tau, (eta, 1.0), colors='k', alpha=0.15)
    # ax.contour(vv, nn, tau, (eta/3, eta, 3*eta), colors='r')
    # ax.contour(vv, nn, tau, (1.0, 3.0), colors='m')
    cs = ax.contour(vv, nn, R0, R0s, linewidths=lws, colors='k')
    clevs = [level for level in clevs_to_label if level in cs.levels]
    ax.clabel(cs, clevs,
              fontsize='x-small', fmt=cformats,
              inline=True, inline_spacing=2, use_clabeltext=True)
    ax.text(62.0, 1e-2, Mlabel, zorder=100, fontsize='x-small', bbox=box_params)
    ax.text(18.0, d["RBW y"][M], RBW_label, rotation=15, fontsize='xx-small', bbox={**box_params, **dict(fc='0.85', ec='0.6')})
    ax.text(16.0, 2e6, RBS_label, rotation=15, fontsize='xx-small', bbox=box_params)
    ax.text(20.0, 300.0, WBS_label, rotation=15, fontsize='xx-small', bbox=box_params)


    #
    # Now do the cooling length
    #
    # Sound speed:
    cs = 11.4*np.sqrt(T0/1e4)
    # pre-shock Mach number
    M0 = vv/cs
    # post-shock Mach number
    M1 = np.sqrt((M0**2 + 3)/(5*M0**2 - 1))
    # post-shock temperature in units of T0
    T1 = (5*M0**2 - 1)*(1 + 3/M0**2) / 16
    # post-shock density
    n1 = nn*4*M0**2 / (M0**2 + 3)
    # post-shock velocity
    v1 = vv*nn/n1
    # Cooling rate
    Lam1 = 3.3e-24 * (T1*T0/1e4)**2.3
    Lam2 = 1e-20 / (T1*T0/1e4)

    k = 3
    Lambda = (Lam1**(-k) + Lam2**(-k))**(-1/k)

    # Heating rate
    Gamma = 1e-26

    # Cooling length in parsec
    dcool = 3*(1e5*v1)*(1.3806503e-16 * T1*T0) / (n1*(Lambda - Gamma)) / 3.085677582e18

    dcool[vv < cs] = np.nan

    # Ratio with respect to adiabatic shell thickness
    h1 = 0.177*R0
    cool_ratio1 = dcool / R0
    # Ratio with respect to isothermal shell thickness
    h2 = 3*R0/(4*M0**2) * (2 / (1 + np.sqrt(1 + (18/M0**2)) ))
    cool_ratio2 = dcool / h2

    cs = ax.contour(vv, nn, cool_ratio1, (1.0,),
                    linewidths=2, colors='b', alpha=0.5)
    ax.clabel(cs, 
              fontsize='xx-small', fmt=r"$d_\mathrm{cool} = R_0$",
              inline=True, inline_spacing=2, use_clabeltext=True)
    cs = ax.contour(vv, nn, cool_ratio2, (1.0,),
                    linewidths=1, colors='b', alpha=0.5)
    ax.clabel(cs, 
              fontsize='xx-small', fmt=r"$d_\mathrm{cool} = h_0$",
              inline=True, inline_spacing=2, use_clabeltext=True)


    ax.set(yscale='log')

axes.set(xlabel=r"$v$, km s$^{-1}$", ylabel=r"$n$, cm$^{-3}$")
sns.despine()
fig.tight_layout()
fig.savefig(figname)

print(figname, end='')
