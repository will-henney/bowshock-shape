import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from bow_projection import (xyprime_t, theta_infinity, theta_0_90,
                            characteristic_radii_projected)
from dragoid_shape import Dragoid

figfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('ticks')
fig, axes = plt.subplots(2, 2, figsize=(8, 8))

# For some reason, exactly 30.0 had problems with R0p
inclinations = [0.0, 15.0, 30.01, 45.0, 60.0, 75.01]
linewidths = [2.4, 2.0, 1.6, 1.2, 0.8, 0.4]
colors = sns.color_palette('magma_r', n_colors=len(inclinations))

# alphas = [1.0, 1.0, 2.0, 2.0]
# mus = [0.05, 0.2, 0.05, 0.2]
alphas = [1.0, 1.0, 4.0, 4.0]
mus = [0.05, 0.2, 0.2, 0.8]
for alpha, mu, ax in zip(alphas, mus, axes.flat):
    shape = Dragoid(alpha=alpha, mu=mu)
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
                label=fr"$i = {inc_dg:.0f}^\circ$",
                color=color, lw=1.5*lw)

    ax.plot([0], [0], 'o', color='k')

    ax.legend(title="Dragoid " + shape.label, ncol=2, loc="center left")
    ax.set(
        xlabel=r"$x' / R_0'$",
        ylabel=r"$y' / R_0'$",
        xlim=[-7, 3],
        ylim=[-5, 5],
    )
    ax.set_aspect('equal', adjustable='box')

sns.despine()
fig.tight_layout()
fig.savefig(figfile)
print(figfile, end='')
