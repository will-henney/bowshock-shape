import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from bow_projection import (xyprime_t, Spline_R_theta_from_function,
                            paraboloid_R_theta, wilkinoid_R_theta,
                            cantoid_R_theta)

figfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('ticks')
fig, axes = plt.subplots(2, 2, figsize=(8, 8))

th = np.linspace(0.0, np.pi, 1001)
inclinations = [0, 15, 30, 45, 60, 75]
linewidths = [2.4, 2.0, 1.6, 1.2, 0.8, 0.4]
colors = sns.color_palette(n_colors=len(inclinations))

for shape_name, ax, R_theta in [
        ["Paraboloid", axes[0, 0],
         Spline_R_theta_from_function(ngrid=301,
                                      shape_func=paraboloid_R_theta)],
        ["Wilkinoid", axes[0, 1],
         Spline_R_theta_from_function(ngrid=101, epsilon=1e-6,
                                      shape_func=wilkinoid_R_theta)],
        [r"Cantoid $\beta = 0.001$", axes[1, 0],
         Spline_R_theta_from_function(ngrid=101, epsilon=1e-6,
                                      shape_func=cantoid_R_theta,
                                      shape_func_pars=(0.001,))],
        [r"Cantoid $\beta = 0.01$", axes[1, 1],
         Spline_R_theta_from_function(ngrid=101, epsilon=1e-6,
                                      shape_func=cantoid_R_theta,
                                      shape_func_pars=(0.01,))],
]:
    for inc_dg, color, lw in zip(inclinations, colors, linewidths):
        inc = np.radians(inc_dg)
        xp, yp = xyprime_t(th, inc, R_theta)
        m = np.isfinite(xp) & np.isfinite(yp)
        if m.sum() == 0:
            # Case of no tangent line at all at this inclination
            continue
        xxp = np.concatenate((xp[m][::-1], xp[m]))
        yyp = np.concatenate((-yp[m][::-1], yp[m]))
        R0p = xxp.max()
        ax.plot(xxp/R0p, yyp/R0p,
                label=fr"$i = {inc_dg:d}^\circ$",
                color=color, lw=lw)

    ax.plot([0], [0], 'o', color='k')

    ax.legend(title=shape_name, ncol=2, loc="center left")
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
