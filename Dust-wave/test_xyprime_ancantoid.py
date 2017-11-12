import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from bow_projection import (xyprime_t, theta_infinity, theta_0_90,
                            characteristic_radii_projected)
from ancantoid_shape import Ancantoid

figfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('ticks')
fig, axes = plt.subplots(2, 2, figsize=(6, 6), sharex=True, sharey=True)

inclinations = [0, 15, 30, 45, 60, 75]
linewidths = [2.4, 2.0, 1.6, 1.2, 0.8, 0.4]
colors = sns.color_palette('magma_r', n_colors=len(inclinations))

for xi, beta, ax in [[0.8, 0.001, axes[0, 0]],
                     [0.8, 0.1, axes[0, 1]],
                     [0.4, 0.001, axes[1, 0]],
                     [0.4, 0.1, axes[1, 1]],]:

    label = "Ancantoid\n" fr"$\beta = {beta:.3f}$, $k = {2/xi - 2:.1f}$"
    shape = Ancantoid(xi=xi, beta=beta)
    th_inf = theta_infinity(shape)
    for inc_dg, color, lw in zip(inclinations, colors, linewidths):
        inc = np.radians(inc_dg)
        th0, th90 = theta_0_90(inc, shape)
        th = np.linspace(th0, th_inf, 301)
        xp, yp = xyprime_t(th, inc, shape)
        m = np.isfinite(xp) & np.isfinite(yp)
        if m.sum() == 0:
            # Case of no tangent line at all at this inclination
            continue
        xxp = np.concatenate((xp[m][::-1], xp[m]))
        yyp = np.concatenate((-yp[m][::-1], yp[m]))
        radii = characteristic_radii_projected(inc, shape)
        R0p = radii['R_0 prime']
        ax.plot(xxp/R0p, yyp/R0p,
                label=fr"$i = {inc_dg:d}^\circ$",
                color=color, lw=1.5*lw)

    ax.plot([0], [0], 'o', color='k')

    ax.legend(title=label, ncol=2, fontsize='small',
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
