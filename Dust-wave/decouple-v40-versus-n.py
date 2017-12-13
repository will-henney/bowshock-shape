import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from vec_root import chandrupatla

figname = sys.argv[0].replace('.py', '.pdf')

sns.set_style('ticks')
sns.set_color_codes('dark')
fig, ax = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(4, 3.5))

# Velocities in units of km/s (10 km/s -> 100 km/s)
vv = 40.0
# Densities in units of 1 pcc  (0.01 -> 1e5)
logngrid = np.linspace(-4.3, 7.3, 600)
nn = 10**logngrid

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

M, L4, eta = 10.0, 0.63, 0.0026

Mlabel = "\n".join([rf"$M = {M:.0f}\, M_\odot$   $\eta = {eta}$",
                    rf"$L = {1e4*L4:.1e}\, L_\odot$".replace("e+0", r"\times 10^")])
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
dv_v = 0.82*(1.0 + alpha)*Zd
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
ax.plot(nn[~m], R0[~m], color="k", lw=0.1)
R0[~m] = np.nan

m = tau <= eta
ax.plot(nn[m], R0[m], label="Wind bow shock")
ax.plot(nn, Rdw, ls="--", label="Dust wave")
ax.plot(nn, Ribs, ls=":", label="Dust-free bow shock")
m = (tau > eta) & (tau < 1.0)
ax.plot(nn[m], R0[m], ls="-.", label="Bow wave")
m = tau >= 1.0
ax.plot(nn[m], R0[m], label="Radiation bow shock")

m = nn < 10.0
ax.plot(nn[m], Rstarstar[m], color='k', lw=0.1)
ax.plot(nn[~m], Rstarstar[~m]*Zd, color='k', lw=0.1)
# ax.contourf(vv, nn, tau, (eta, 1.0), colors='k', alpha=0.2)
# ax.contourf(vv, nn, (~m)*np.ones_like(m).astype(float), (0.5, 1.5), colors='c', alpha=0.2)
# if np.any(~m):
#     alpha[m] = np.nan
#     cs = ax.contour(vv, nn, alpha, (0.25, 0.5, 1.0, 2.0, 4.0, 20.0, 50.0),
#                     linewidths=0.3,
#                     colors='m', alpha=0.5)
#     ax.clabel(cs,
#               fontsize='x-small', fmt=r"$\alpha = %.2f$",
#               inline=True, inline_spacing=0.2, use_clabeltext=True)
# cs = ax.contour(vv, nn, R0, R0s, linewidths=lws, colors='k')
# clevs = [level for level in clevs_to_label if level in cs.levels]
# ax.clabel(cs, clevs,
#           fontsize='xx-small', fmt=cformats,
#           inline=True, inline_spacing=0.2, use_clabeltext=True)


# cs = ax.contour(vv, nn, Rdw, R0s, linewidths=lws, linestyles='dashed', colors='y')
# clevs = [level for level in clevs_to_label if level in cs.levels]
# ax.clabel(cs, clevs,
#           fontsize='xx-small', fmt=cformats, colors='y',
#           inline=True, inline_spacing=0.2, use_clabeltext=True)
# cs = ax.contour(vv, nn, Ribs, R0s, linewidths=lws, colors='k', alpha=0.5)
# clevs = [level for level in clevs_to_label if level in cs.levels]
# ax.clabel(cs, clevs,
#           fontsize='xx-small', fmt=cformats, alpha=0.5, 
#           inline=True, inline_spacing=0.2, use_clabeltext=True)

ax.text(3e-5, 9e-5, "Stellar parameters:\n" + Mlabel,
        linespacing=1.5, fontsize='small', bbox=box_params)
ax.text(1e6, 9e-3, "Relative speed:\n$v = 40$ km/s", ha="right",
        linespacing=1.5, fontsize='small', bbox=box_params)
ax.text(5.0, 0.15, "$R_{_{**}}$",
        fontsize='small', ha="center", va="center", bbox=box_params)
ax.text(1e6, 0.15/100, "$Z_\mathrm{d}\,R_{_{**}} $",
        fontsize='small', ha="center", va="center", bbox=box_params)

ax.text(3e-5, 1.2e-3,
        "Dust parameters:\n" + f"$Z_\mathrm{{d}} = {Zd}$"
        + fr"     $\xi = {xi}$" + "\n"
        + f"$\kappa_\mathrm{{d}} = {kappa}$ cm$^2$/g",
        linespacing=1.5, fontsize='small', bbox=box_params)


leg = ax.legend(loc="upper right", frameon=True, fontsize="small")

ax.set(xscale='log', yscale='log')

ax.set(
    ylabel=r"Bow radius: $R$, pc",
    xlabel=r"Ambient density: $n$, cm$^{-3}$",
    #ylim=[1.5e-4, 3.0],
)
sns.despine()
fig.tight_layout()
fig.savefig(figname)

print(figname, end='')
