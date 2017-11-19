import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import bow_projection as bp
from dragoid_shape import Dragoid

figfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('ticks')
fig, ax = plt.subplots(figsize=(10, 4))

th = np.linspace(-np.pi, np.pi, 1001)
th_dg = np.degrees(th)

alpha_drags = [0.25, 1.0] + [1.0, 1.0, 4.0, 4.0]
mus = [None]*2 + [.05, 0.2, 0.2, 0.8]

for alpha_drag, mu in zip(alpha_drags, mus):
    shape = Dragoid(alpha=alpha_drag, mu=mu, lowess_frac=0.1)
    alpha_slope = np.degrees(bp.alpha(shape.thgrid, shape))
    fac = (90.0 - alpha_slope)/np.degrees(shape.thgrid)
    ax.plot(np.degrees(shape.thgrid), fac,
            lw=0.8, label=shape.label)

ax.legend(title=r"Dragoid shapes")
ax.set(
    xlabel=r"Polar angle: $\theta$, degrees",
    ylabel=r"Slope angle, $(90 - \alpha)/\theta$",
    xlim=[0, 180],
    ylim=[0, 1.1],
    xticks=[0, 30, 60, 90, 120, 150, 180],
)
sns.despine()
fig.tight_layout()
fig.savefig(figfile)
print(figfile, end='')
