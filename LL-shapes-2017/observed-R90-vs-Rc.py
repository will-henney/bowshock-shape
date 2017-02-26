import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker
import seaborn as sns
import json

plotfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('white')

fig, ax = plt.subplots(figsize=(5, 5))
Rc_grid = np.linspace(0.0, 10.0, 2000)
R90_T0_grid = np.sqrt(2*Rc_grid)
R90_T1_grid = np.sqrt(2*Rc_grid - 1.0)
R90_T1_grid[~np.isfinite(R90_T1_grid)] = 0.0 

ax.fill_between(Rc_grid, R90_T1_grid, R90_T0_grid, color='k', alpha=0.2)
ax.fill_between(Rc_grid, R90_T0_grid, color='k', alpha=0.1)
ax.plot(Rc_grid, R90_T0_grid, c='k', lw=0.5)
ax.axhline(1.0, lw=0.5, alpha=0.5, color='k', zorder=-1)
ax.axvline(1.0, lw=0.5, alpha=0.5, color='k', zorder=-1)
ax.plot([0.0, 10.0], [0.0, 10.0], lw=0.5, alpha=0.5, color='k', zorder=-1)

# Plot the observational data
arc_data_ll = '../read-shapes-LL/radii-set.json'
data = json.load(open(arc_data_ll))['outer']
sources = list(data['R0'].keys())
x, y = [], []
for source in sources:
    try:
        y.append(0.5*(data['R90'][source] + data['Rm90'][source])
                 / data['R0'][source])
    except:
        try:
            y.append(data['R90'][source] / data['R0'][source])
        except:
            try:
                y.append(data['Rm90'][source] / data['R0'][source])
            except:
                continue
    x.append(data['Rc'][source] / data['R0'][source])

ax.scatter(x, y)

ax.set(
    yscale='linear',
    xscale='linear',
    xlim=[0.0, 5.1],
    ylim=[0.0, 5.1],
#    ylim=[-3.0, 1.1],
    xlabel=r"Projected dimensionless radius of curvature: $\widetilde{R}_{c}{}'$",
    ylabel=r"Projected dimensionless perpendicular radius: $\widetilde{R}_{90}{}'$",
)        


fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
