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

Simulation.lowess_frac = 0.1
mode = "positive"
for label, shape in [
        ["MHD open",
         Simulation(name="M17-MHD2040-AllB7", force_open=True, mode=mode)],
        ["MHD closed",
         Simulation(name="M17-MHD2040-AllB7", force_open=False, mode=mode)],
        ["HD open",
         Simulation(name="M17-HD2040", force_open=True, mode=mode)],
        ["HD closed",
         Simulation(name="M17-HD2040", force_open=False, mode=mode)],
]:
    ax.plot(th_dg, omega(th, shape), label=label)

ax.legend(title=r"Simulation shapes")
ax.axhline(1.0, xmin=0.35, xmax=0.65, color='white', lw=4, zorder=100)
ax.axhline(1.0, xmin=0.35, xmax=0.65, color='k', lw=1, ls=':', zorder=101)
ax.axhspan(0.0, 1.0, color='k', alpha=0.05, ec='none')
ax.set_yscale('symlog', linthreshy=1.0, linscaley=0.5)
ax.set(
    xlabel=r"Polar angle: $\theta$, degrees",
    ylabel=r"$\omega \equiv R^{-1} d R / d \theta$",
    xlim=[0, 180],
    ylim=[-0.5, 10.1],
    xticks=[0, 30, 60, 90, 120, 150, 180],
)
sns.despine()
fig.tight_layout()
fig.savefig(figfile)
print(figfile, end='')
