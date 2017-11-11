import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker
import seaborn as sns
import conic_parameters
from equation6 import Shell

plotfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('ticks')
sns.set_color_codes(palette='deep')
fig, ax = plt.subplots(figsize=(4, 5))

ks = [None, 0.0, 0.5, 3.0, 8.0]
colors = 'kkrmb'
linestyles = ['--'] + ['-']*4
betas = [0.01, 0.1, 0.5]
lws = [2.0, 1.5, 1.0]
theta = np.linspace(0,0.99*np.pi)

for k, color, ls in zip(ks, colors, linestyles):
    for beta, lw in zip(betas, lws):
        if k is None:
            shell = Shell(beta=beta, innertype="isotropic")
            label = fr"Cantoid"
        else:
            xi = 2.0 / (k + 2.0)
            shell = Shell(beta=beta, innertype="anisotropic", xi=xi)
            label = fr"Ancantoid, $k = {k:.1f}$"
        if lw < 2.0:
            label = "_nolabel_"

        R = shell.radius(theta)
        R[R < 0.0] = np.nan

        ax.plot(R*np.cos(theta), R*np.sin(theta),
                ls=ls, color=color, lw=lw, label=label)
        ax.plot(R*np.cos(theta), -R*np.sin(theta),
                ls=ls, color=color, lw=lw, label="_nolabel_")

ax.axhline(0.0, color='k', ls=':')
ax.axvline(0.0, color='k', ls=':')
ax.plot([0.0, 1.0], [0.0, 0.0], 'o', color='k')

for xtext, beta in zip([-0.2, 0.13, 0.4], betas):
    ax.text(xtext, -0.2, fr"$\beta = {beta:.2f}$", fontsize='small',
            ha='center', bbox={'fc': 'white', 'ec': 'none', 'pad': 2})

ax.legend(ncol=1, fontsize='small', frameon=True,
          labelspacing=0.8, loc='upper right', title=None)
ax.set(
    xlim=[-0.45, 1.05],
    ylim=[-0.25, 1.85],
    xlabel=r"$z / D$",
    ylabel=r"$r / D$",
    xticks=[-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
    yticks=[-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8],
)
ax.set_aspect('equal')
sns.despine(trim=True)
fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
