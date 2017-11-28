import sys
import numpy as np
from simulation_shape import Simulation
sys.path.append("../Dust-wave")
from bow_projection import omega
from matplotlib import pyplot as plt
import seaborn as sns

figfile = sys.argv[0].replace('.py', '.pdf')


sns.set_style('ticks')
fig, ax = plt.subplots()

th = np.linspace(-np.pi, np.pi, 1001)
th_dg = np.degrees(th)

for sim in ["M17-MHD2040-AllB7", "M17-HD2040"]:
    shape = Simulation(name=sim)
    ax.plot(th_dg, omega(th, shape), label=sim)

ax.legend(title=r"Simulation shapes")
ax.axhline(1.0, xmin=0.35, xmax=0.65, color='white', lw=4, zorder=100)
ax.axhline(1.0, xmin=0.35, xmax=0.65, color='k', lw=1, ls=':', zorder=101)
ax.axhspan(0.0, 1.0, color='k', alpha=0.05, ec='none')
ax.set_yscale('symlog', linthreshy=1.0, linscaley=0.5)
ax.set(
    xlabel=r"Polar angle: $\theta$, degrees",
    ylabel=r"$\omega \equiv R^{-1} d R / d \theta$",
    xlim=[0, 180],
    ylim=[0.0, 2e2],
    xticks=[0, 30, 60, 90, 120, 150, 180],
)
sns.despine()
fig.tight_layout()
fig.savefig(figfile)
print(figfile, end='')
