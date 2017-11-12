import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from bow_projection import (xyprime_t, theta_infinity, theta_0_90,
                            paraboloid_R_theta, wilkinoid_R_theta,
                            cantoid_R_theta)

figfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('ticks')
fig, axes = plt.subplots(2, 2, figsize=(6, 6), sharex=True, sharey=True)

inclinations = [0, 15, 30, 45, 60.1, 75]
linewidths = [2.4, 2.0, 1.6, 1.2, 0.8, 0.4]
colors = sns.color_palette('magma_r', n_colors=len(inclinations))

for shape_name, ax, R_theta, extra_pars in [
        ["Paraboloid", axes[0, 0], paraboloid_R_theta, ()],
        ["Wilkinoid", axes[0, 1], wilkinoid_R_theta, ()],
        ["Cantoid\n" r"$\beta = 0.001$", axes[1, 0], cantoid_R_theta, (0.001,)],
        ["Cantoid\n" r"$\beta = 0.01$", axes[1, 1], cantoid_R_theta, (0.01,)],
]:
    th_inf = theta_infinity(R_theta, *extra_pars)
    for inc_dg, color, lw in zip(inclinations, colors, linewidths):
        inc = np.radians(inc_dg)
        th0, th90 = theta_0_90(inc, R_theta, *extra_pars)
        th = np.linspace(th0, th_inf, 101)
        xp, yp = xyprime_t(th, inc, R_theta, *extra_pars)
        m = np.isfinite(xp) & np.isfinite(yp)
        if m.sum() == 0:
            # Case of no tangent line at all at this inclination
            continue
        xxp = np.concatenate((xp[m][::-1], xp[m]))
        yyp = np.concatenate((-yp[m][::-1], yp[m]))
        R0p = xxp.max()
        ax.plot(xxp/R0p, yyp/R0p, label=fr"$i = {int(inc_dg):d}^\circ$", lw=1.5*lw, color=color)

    ax.plot([0], [0], 'o', color='k')

    ax.legend(title=shape_name, ncol=2, fontsize='small',
              handlelength=1.0, handletextpad=0.5, columnspacing=0.3,
              loc="center left")
    ax.set_aspect('equal', adjustable='box-forced')

axes[-1,0].set(
    xlabel=r"$x' / R_0'$",
    ylabel=r"$y' / R_0'$",
    xlim=[-7, 3],
    ylim=[-5, 5],
)
sns.despine()
fig.tight_layout(pad=0.3, h_pad=0.1, w_pad=0.1)
fig.savefig(figfile)
print(figfile, end='')
