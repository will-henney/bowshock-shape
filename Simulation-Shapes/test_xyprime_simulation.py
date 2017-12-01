import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
sys.path.append("../Dust-wave")
from bow_projection import (xyprime_t, theta_infinity, theta_0_90,
                            characteristic_radii_projected)
from simulation_shape import Simulation

figfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('ticks')
fig, axes = plt.subplots(2, 2, figsize=(6, 6), sharex=True, sharey=True)

# For some reason, exactly 30.0 had problems with R0p
inclinations = [0.0, 15.0, 30.01, 45.0, 60.0, 75.01]
linewidths = [2.4, 2.0, 1.6, 1.2, 0.8, 0.4]
colors = sns.color_palette('magma_r', n_colors=len(inclinations))

sims = ["M17-MHD2040-AllB7", "M17-HD2040", "M17-MHD2040-AllB7", "M17-HD2040"]

labels = ["MHD open", "MHD closed", "HD open", "HD closed"]
mode = "negative"
shapes = [
    Simulation(name="M17-MHD2040-AllB7", force_open=True, cheby_degree=12),
    Simulation(name="M17-MHD2040-AllB7", force_open=False, cheby_degree=12),
    Simulation(name="M17-HD2040", force_open=True, cheby_degree=12, extrap_degree=1),
    Simulation(name="M17-HD2040", force_open=False, cheby_degree=12, extrap_degree=1)]

for label, shape, ax in zip(labels, shapes, axes.flat):
    th_inf = theta_infinity(shape)
    th_inf = max(th_inf, np.pi)
    for inc_dg, color, lw in zip(inclinations, colors, linewidths):
        inc = np.radians(inc_dg)
        th0, th90 = theta_0_90(inc, shape)
        if not np.isfinite(th0):
            th0 = 0.0
        th = np.linspace(th0, th_inf, 301)
        xp, yp = xyprime_t(th, inc, shape)
        m = np.isfinite(xp) & np.isfinite(yp)
        # if m.sum() == 0:
        #     # Case of no tangent line at all at this inclination
        #     continue
        xxp = np.concatenate((xp[m][::-1], xp[m]))
        yyp = np.concatenate((-yp[m][::-1], yp[m]))
        radii = characteristic_radii_projected(inc, shape)        
        R0p = radii['R_0 prime']
        ax.plot(xxp/R0p, yyp/R0p,
                label=fr"$i = {inc_dg:.0f}^\circ$",
                color=color, lw=1.5*lw)

    ax.plot([0], [0], 'o', color='k')

    if "open" in label:
        ax.legend(fontsize='small',
                  handlelength=1.0, handletextpad=0.5, columnspacing=0.3,
                  ncol=2, loc="center left")
    ax.text(0, 4, label, ha='right', va='top')
    ax.set_aspect('equal', adjustable='box-forced')

axes[-1,0].set(
    xlabel=r"$x' / R_0'$",
    ylabel=r"$y' / R_0'$",
    xlim=[-7, 3],
    ylim=[-5, 5],
)

sns.despine()
for ax in axes.flat:
    ax.label_outer()
fig.tight_layout(pad=0.3, h_pad=0.1, w_pad=0.1)
fig.savefig(figfile)
print(figfile, end='')
