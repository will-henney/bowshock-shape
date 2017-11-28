import sys
import numpy as np
from simulation_shape import Simulation
from matplotlib import pyplot as plt
import seaborn as sns

figfile = sys.argv[0].replace('.py', '.pdf')


sns.set_style('ticks')
fig, ax = plt.subplots()

th = np.linspace(-np.pi, np.pi, 1001)
th_dg = np.degrees(th)

for sim in ["M17-MHD2040-AllB7", "M17-HD2040"]:
    shape = Simulation(name=sim)
    ax.plot(np.degrees(shape.thgrid), shape.Rgrid,
            color='b', alpha=0.2, lw=2, label='_nolabel_')
    ax.plot(th_dg, shape(th), lw=0.8, label=sim)

ax.legend(title=r"Simulation shapes")
ax.set(
    xlabel=r"Polar angle: $\theta$, degrees",
    ylabel=r"$R$",
    xlim=[0, 180],
    yscale='log',
    ylim=[0.9, 200.0],
    xticks=[0, 30, 60, 90, 120, 150, 180],
)
sns.despine()
fig.tight_layout()
fig.savefig(figfile)
print(figfile, end='')
