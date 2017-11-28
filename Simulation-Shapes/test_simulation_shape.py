import sys
import numpy as np
from simulation_shape import Simulation
from matplotlib import pyplot as plt
import seaborn as sns

figfile = sys.argv[0].replace('.py', '.pdf')


sns.set_style('ticks')
fig, ax = plt.subplots(figsize=(4, 4))

th = np.linspace(-np.pi, np.pi, 1001)
th_dg = np.degrees(th)

for label, shape in [
        ["MHD open",
         Simulation(name="M17-MHD2040-AllB7",
                    force_open=True, extrap_degree=2)],
        ["MHD closed",
         Simulation(name="M17-MHD2040-AllB7",
                    force_open=False, extrap_degree=2)],
        ["HD open",
         Simulation(name="M17-HD2040",
                    force_open=True, extrap_degree=1)],
        ["HD closed",
         Simulation(name="M17-HD2040",
                    force_open=False, extrap_degree=1)],
]:
    ax.plot(np.degrees(shape.thgrid), shape.Rgrid,
            color='b', alpha=0.2, lw=2, label='_nolabel_')
    ax.plot(th_dg, shape(th), lw=0.8, label=label)

ax.legend(title=r"Simulation shapes")
ax.set(
    xlabel=r"Polar angle: $\theta$, degrees",
    ylabel=r"$R$",
    xlim=[-180, 180],
    yscale='log',
    ylim=[0.9, 30.0],
    xticks=[0, 30, 60, 90, 120, 150, 180],
)
sns.despine()
fig.tight_layout()
fig.savefig(figfile)
print(figfile, end='')
