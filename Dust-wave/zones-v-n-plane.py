import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from vec_root import chandrupatla

figname = sys.argv[0].replace('.py', '.pdf')

sns.set_style('ticks')
sns.set_color_codes('dark')
fig, axes = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(8, 4))
stardata = [
    [10.0, 0.63, 0.0066, 1e-4, axes[0]],
    [20.0, 5.453, 0.1199, 0.1, axes[1]],
    [40.0, 22.19, 0.4468, 1.0, axes[2]],
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
    "RBW y": {10.0: 4000.0, 20.0: 5000.0, 40.0: 2500.0},
    "trapped y": {10.0: 3.5e4, 20.0: 2.5e5, 40.0: 1.5e5},
    "trapped bg": {10.0: '0.85', 20.0: 'w', 40.0: 'w'},
    "IF tau": {10.0: 0.2, 20.0: 3.7, 40.0: 6.4},
    "IF tau gas": {10.0: 5.0, 20.0: 5.0, 40.0: 5.0},
}

T0 = 8000
kappa = 600.0
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

    # Ionization parameter - fiducial
    U = 2.789*S49 / (R0**2 * nn)
    # Shell ionization parameter - assume compression by M^2
    Ush = U/(vv/10)**2
    # Ionization fraction in shell
    ysh = 1.0 - 1.0/(3.5e5*Ush)

    # cs = ax.contour(vv, nn, ysh,
    #                 (0.9, 0.99, 0.999, 0.9999, 0.99999),
    #                 linewidths=1.0,
    #                 colors='g', alpha=0.5)
    # ax.clabel(cs,
    #           fontsize=5, colors='g', fmt='$%.5f$', 
    #           inline=True, inline_spacing=1, use_clabeltext=True)
    # ax.text(75, 6e-3, fr"$y_\mathrm{{in}} = {ysh.max():.5f}$",
    #         fontsize=5, color='g', alpha=0.5)
    # ax.text(12, 1.3e6, fr"$y_\mathrm{{in}} = {ysh.min():.5f}$",
    #         fontsize=5, color='g', alpha=0.5)

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

    # y^2 / (1 - y) = CU
    # y = 0.1 => y^2 / (1 - y) = 0.0111111111111
    # y = 0.5 => y^2 / (1 - y) = 0.5
    # y = 0.9 => y^2 / (1 - y) = 8.1
    # y = 0.99 => y^2 / (1 - y) = 98.1

    # cu = 3.5e5*Uout

    # y_out = 0.5*cu * (np.sqrt(1.0 + 4/cu) - 1.0)
    # cs = ax.contour(vv, nn, Uout,
    #                 np.logspace(-7.0, 1.0, 9),
    #                 linewidths=np.linspace(0.2, 1.5, 9),
    #                 colors='r', alpha=0.5)
    # ax.clabel(cs,x
    #           fontsize='xx-small', colors='r', fmt='%.0e', 
    #           inline=True, inline_spacing=1, use_clabeltext=True)

    #ax.contour(vv, nn, 3.5e5*Uout, (0.5, 98.1), colors='r', alpha=0.5)
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

    # ax.contour(vv, nn, tau, tau_dust_IF,
    #            colors='k', linestyles='--', linewidths=0.8)
    # ax.contour(vv, nn, tau_gas, tau_gas_IF,
    #            colors='y', linewidths=0.4)
    arrows = r"$\uparrow\!\!\!\!\uparrow$"
    ax.text(60, d["trapped y"][M],
            rf"{arrows} Trapped i-front, $\tau_\mathrm{{d}} = {tau_dust_IF:.1f}$, $\tau_\mathrm{{gas}} = {tau_gas_IF:.1f}$ {arrows}",
            ha='center', va='center',
            fontsize='xx-small', color='r', alpha=0.5, rotation=10,
            bbox=dict(fc=d["trapped bg"][M], ec='none', pad=0.1)
    )
    # ax.clabel(cs, (0.5,),
    #           manual=((60, 1e4),),
    #           fontsize='xx-small', colors='r',
    #           fmt=r'$\uparrow\uparrow$ Trapped i-front $\uparrow\uparrow$', 
    #           inline=True, inline_spacing=10, use_clabeltext=True)

    ax.contourf(vv, nn, tau, (eta, 1.0), colors='k', alpha=0.15)
    # ax.contour(vv, nn, tau, (eta/3, eta, 3*eta), colors='r')
    # ax.contour(vv, nn, tau, (1.0, 3.0), colors='m')
    cs = ax.contour(vv, nn, R0, R0s, linewidths=lws, colors='k')
    clevs = [level for level in clevs_to_label if level in cs.levels]
    ax.clabel(cs, clevs,
              fontsize='x-small', fmt=cformats,
              inline=True, inline_spacing=2, use_clabeltext=True)
    ax.text(18.0, 3e-3, Mlabel, zorder=100, fontsize='x-small', bbox=box_params)
    ax.text(18.0, d["RBW y"][M], RBW_label, rotation=15, fontsize='xx-small', bbox={**box_params, **dict(fc='0.85', ec='0.6')})
    ax.text(16.0, 1e6, RBS_label, rotation=15, fontsize='xx-small', bbox=box_params)
    ax.text(20.0, 15.0, WBS_label, rotation=15, fontsize='xx-small', bbox=box_params)


    #
    # Now do the cooling length
    #
    # pre-shock Mach number
    M0 = vv/10.0
    # post-shock Mach number
    M1 = np.sqrt((M0**2 + 3)/(5*M0**2 - 1))
    # post-shock temperature in units of T0
    T1 = (5*M0**2 - 1)*(1 + 3/M0**2) / 16
    # post-shock density
    n1 = nn*4*M0**2 / (M0**2 + 3)
    # post-shock velocity
    v1 = vv*nn/n1
    # Cooling rate
    Lam1 = 3.3e-24 * T1**2.3
    Lam2 = 1e-20 / T1
    k = 3
    Lambda = (Lam1**(-k) + Lam2**(-k))**(-1/k)
    # Cooling length in parsec
    dcool = 3*(1e5*v1)*(1.3806503e-16 * T1*T0) / (n1*Lambda) / 3.085677582e18

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

axes[0].set(xlabel=r"$v$, km s$^{-1}$", ylabel=r"$n$, cm$^{-3}$")
sns.despine()
for ax in axes:
    ax.label_outer()
fig.tight_layout()
fig.savefig(figname)

print(figname, end='')
