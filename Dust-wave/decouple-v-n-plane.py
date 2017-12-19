import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from vec_root import chandrupatla

figname = sys.argv[0].replace('.py', '.pdf')

sns.set_style('ticks')
sns.set_color_codes('dark')
fig, ax = plt.subplots(1, 1, figsize=(3, 4))
stardata = [
    [10.0, 0.63, 0.0066, ax],
]

# Velocities in units of km/s (10 km/s -> 100 km/s)
vgrid = np.linspace(10.0, 100.0, 600)
# Densities in units of 1 pcc  (0.01 -> 1e5)
logngrid = np.linspace(-3.3, 6.3, 600)
# 2d versions of velocity and density grids
vv, nn = np.meshgrid(vgrid, 10**logngrid)

def rstar(v10, n, L4, c10=0.0):
    """Characteristic radius in pc"""
    return 2.21*np.sqrt(L4/(n*(v10**2 + c10**2)))

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

# Dust-gas mass ratio
Zd = 0.01
kappa = 600.0
# Decoupling efficiency: 2 Q_p / Q_drag
xi = 0.07

for M, L4, eta, ax in stardata:
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
    # Drag coefficient
    alpha = (2.0/xi**0.5) * ts / Zd

    # Dust wave radius
    Rdw = Rstarstar / (1.0 + alpha)
    # Rdw = xi**0.5 * Rs

    # Pure wind bow shock radius
    Rwbs = eta**0.5 * Rs

    # Gas deceleration through dust wave Delta v/v
    dv_v = 0.75*(1.0 + alpha)*Zd
    # Inner bow shock radius inside a dust wave
    Ribs = eta**0.5 * rstar((1.0 - dv_v)*vv/10, nn, L4)
    # Don't let it be bigger than dust wave
    m = Ribs > Rdw
    Ribs[m] = Rdw[m]

    # Remove cases where dust wave will not exist
    m = Rdw < R0
    Rdw[m] = np.nan 
    Ribs[m] = np.nan 

    # Choose only the R0 where there is not also a dust wave
    R0[~m] = np.nan

    ax.contourf(vv, nn, tau, (eta, 1.0), colors='k', alpha=0.2)
    ax.contourf(vv, nn, (~m)*np.ones_like(m).astype(float), (0.5, 1.5), colors='c', alpha=0.2)
    if np.any(~m):
        alpha[m] = np.nan
        #ax.contourf(vv, nn, alpha, (0.5, 2.0), colors='m', alpha=0.2)
        cs = ax.contour(vv, nn, alpha, (0.5, 1.0, 2.0, 20.0, 80.0),
                        linewidths=0.3,
                        colors='m', alpha=0.5)
        ax.clabel(cs,
                  fontsize='xx-small', fmt=r"$\alpha = %.2f$",
                  inline=True, inline_spacing=2, use_clabeltext=True)
    # ax.contour(vv, nn, tau, (eta/3, eta, 3*eta), colors='r')
    # ax.contour(vv, nn, tau, (1.0, 3.0), colors='m')
    cs = ax.contour(vv, nn, R0, R0s, linewidths=lws, colors='k')
    clevs = [level for level in clevs_to_label if level in cs.levels]
    ax.clabel(cs, (0.001, 0.01, 0.1, 1.0),
              manual=((25, 5e5), (15, 3e3), (85, 1.2e-2), (20, 2e-3)),
              fontsize='xx-small', fmt=cformats,
              inline=True, inline_spacing=8, use_clabeltext=False)


    # cs = ax.contour(vv, nn, Rstarstar, R0s, linewidths=lws, colors='m', alpha=0.5)
    # clevs = [level for level in clevs_to_label if level in cs.levels]
    # ax.clabel(cs, clevs,
    #           fontsize='small', fmt=cformats,
    #           inline=True, inline_spacing=2, use_clabeltext=True, alpha=0.5, colors='m')
    cs = ax.contour(vv, nn, Rdw, R0s, linewidths=lws, linestyles='dashed', colors='y')
    clevs = [level for level in clevs_to_label if level in cs.levels]
    ax.clabel(cs, (0.001, 0.01, 0.1),
              fontsize='xx-small', fmt=cformats, colors='y',
              inline=True, inline_spacing=0.2, use_clabeltext=True)
    cs = ax.contour(vv, nn, Ribs, R0s, linewidths=lws, colors='k', alpha=0.5)
    ax.clabel(cs, (0.001, 0.01),
              manual=((73, 2e2), (82, 2e0)),
              fontsize='xx-small', fmt=cformats, alpha=0.5, 
              inline=True, inline_spacing=8, use_clabeltext=False)

    ax.text(30.0, 1.5e-3, Mlabel, fontsize='small', bbox=box_params)
    # ax.text(18.0, 4000.0, RBW_label, rotation=15, fontsize='small', bbox={**box_params, **dict(fc='0.8', ec='0.6')})
    # ax.text(18.0, 2e6, RBS_label, rotation=15, fontsize='small', bbox=box_params)
    # ax.text(18.0, 30.0, WBS_label, rotation=15, fontsize='small', bbox=box_params)



    ax.set(yscale='log')

ax.set(xlabel=r"$v$, km s$^{-1}$", ylabel=r"$n$, cm$^{-3}$")
sns.despine()
fig.tight_layout()
fig.savefig(figname)

print(figname, end='')
