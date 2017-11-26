import sys
import json
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from equation6 import Shell

plotfile = sys.argv[0].replace('.py', '.pdf')

def compensate(R, theta):
    """Compensated inversion of R(theta)"""
    return 1.0/R - 0.5*(1 + np.cos(theta))

sns.set_style('ticks')
sns.set_color_codes(palette='deep')
fig, ax = plt.subplots()

# Plot confocal parabola
ax.axhline(0.0, ls='-', c='k', lw=0.5)

# Plot wilkinoid
mugrid = np.linspace(-1.0, 1.0, 100)[::-1]
thgrid = np.arccos(mugrid)
# ax.plot(mugrid, compensate(bp.wilkinoid_R_theta(thgrid), thgrid),
#         '-', c='k', lw=1.5)

betas = [0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 0.0003, 0.0001, 1e-6]

ks = [None, 0.0, 0.5]
colors = '0.5', 'k', 'r'

# Plot cantoids
for k, color in zip(ks, colors):
    for beta in betas:
        label = "_nolabel_"
        if k is None:
            shell = Shell(beta=beta, innertype="isotropic")
            if beta == betas[0]:
                label = fr"Cantoid"
        else:
            xi = 2.0 / (k + 2.0)
            shell = Shell(beta=beta, innertype="anisotropic", xi=xi)
            if beta == betas[0]:
                label = fr"Ancantoid, $k = {k:.1f}$"
        R = shell.radius(thgrid)
        R[R < 0.0] = np.nan
        R /= shell.R0
        ax.plot(mugrid, compensate(R, thgrid), '-', c=color, lw=0.5, label=label)

# Fill in forbidden zone
ax.fill_between(mugrid, -0.5*(1.0 + mugrid), -1.0, color='k', alpha=0.4)

ax.legend()

ax.set(
    xlim=[-1.02, 1.02],
    ylim=[-0.255, 0.255],
    xlabel=r"$\cos \,\theta$",
    ylabel=r"Departure function, $\Delta(\cos\theta)$",
)
sns.despine(trim=True)
fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
