import sys
import json
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from equation6 import Shell
sys.path.append("../Dust-wave")
import bow_projection as bp
import ancantoid_shape

plotfile = sys.argv[0].replace('.py', '.pdf')

def compensate(R, theta):
    """Compensated inversion of R(theta)"""
    return 1.0/R - 0.5*(1 + np.cos(theta))

sns.set_style('ticks')
sns.set_color_codes(palette='deep')
fig, ax = plt.subplots(figsize=(4, 4))

# Plot confocal parabola
ax.axhline(0.0, ls=':', c='k', lw=0.5)

gray = '0.6'

# Plot wilkinoid
mugrid = np.linspace(-1.0, 1.0, 100)[::-1]
thgrid = np.arccos(mugrid)
ax.plot(mugrid, compensate(bp.wilkinoid_R_theta(thgrid), thgrid),
        '-', c=gray, lw=2)

betas = [0.1, 0.01, 0.001, 0.0001, 1e-7]

ks = [None, 0.0, 0.5]
colors = gray, 'k', 'r'

# Plot cantoids
for k, color in zip(ks, colors):
    for beta in betas:
        if k is None:
            R = bp.cantoid_R_theta(thgrid, beta)
            shape_label = fr"Cantoid/Wilkinoid"
        else:
            xi = 2.0 / (k + 2.0)
            shape = ancantoid_shape.Ancantoid(xi=xi, beta=beta, n=301)
            R = shape(thgrid)
            shape_label = fr"Ancantoid $k = {k:.1f}$"
        label = shape_label if beta == betas[-1] else "_nolabel_"
        lw = 1.5 if beta == betas[-1] else 0.5
        ax.plot(mugrid, compensate(R, thgrid),
                '-', c=color, lw=lw, label=label)

# Fill in forbidden zone
#ax.fill_between(mugrid, -0.5*(1.0 + mugrid), -1.0, color='k', alpha=0.4)

ax.annotate(r"$\beta = 0.1$", (-0.6, -0.15), fontsize='small', ha='right')
ax.annotate(r"$\beta = 0.01$", (-0.63, -0.04), fontsize='small', ha='right')
#ax.annotate(r"$\beta = 0.001$", (-0.6, -0.2), fontsize='small', ha='right')

ax.annotate(r"$\beta = 0$", (-0.48, 0.185), color='r', fontsize='small', ha='center')
ax.annotate(r"$\beta = 0$", (-0.56, 0.12), fontsize='small', ha='center')

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
