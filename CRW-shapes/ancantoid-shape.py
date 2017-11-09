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
fig, ax = plt.subplots(figsize=(6, 5))

ks = [None, 0.5, 3.0, 8.0]
colors = 'krmb'
betas = [0.01, 0.1, 0.5]
lws = [2, 1.5, 1]
theta = np.linspace(0,0.99*np.pi)

for k, color in zip(ks, colors):
    for beta, lw in zip(betas, lws):
        if k is None:
            shell = Shell(beta=beta, innertype="isotropic")
            label = fr"Cantoid, $\beta = {beta:.2f}$"
        else:
            xi = 2.0 / (k + 2.0)
            shell = Shell(beta=beta, innertype="anisotropic", xi=xi)
            label = fr"Ancantoid, $k = {k:.1f}$, $\beta = {beta:.2f}$"

        R = shell.radius(theta)
        R[R < 0.0] = np.nan

        ax.plot(R*np.cos(theta), R*np.sin(theta), color=color, lw=lw, label=label)
        ax.plot(R*np.cos(theta), -R*np.sin(theta), color=color, lw=lw, label="_nolabel_")

ax.axhline(0.0, color='k', ls='--')
ax.plot([0.0, 1.0], [0.0, 0.0], 'o', color='k')
ax.legend(ncol=1, fontsize='small', frameon=True,
          labelspacing=0.8, loc='upper right', title='Shape')
ax.set(
    yscale='linear',
    xlim=[-0.45, 1.18],
    ylim=[-0.25, 1.05],
    xlabel=r"$z / D$",
    ylabel=r"$r / D$",
)
ax.set_aspect('equal')
sns.despine(trim=True)
fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
