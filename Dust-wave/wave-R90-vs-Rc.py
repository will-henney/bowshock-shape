import sys
import json
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker
import seaborn as sns

try:
    WAVE_ID = sys.argv[1]
    BASE_SHAPE_IDS = sys.argv[2:]
except:
    sys.exit(f"Usage: {sys.argv[0]} WAVE_ID SHAPE_ID [SHAPE_ID ...]")

plotfile = sys.argv[0].replace('.py', f'-{WAVE_ID}.pdf')

sns.set_style('ticks')
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

# Put a cross at the Wilkinoid coordinates: [5/3, sqrt(3)]
ax.plot([5./3.], [np.sqrt(3.0)], '+', c='w', ms=10, alpha=1.0)

nshapes = len(BASE_SHAPE_IDS)
colors = sns.color_palette(n_colors=nshapes)
for base_shape_id, color in zip(BASE_SHAPE_IDS, colors):
    prefix = f"{base_shape_id}-wave-{WAVE_ID}"
    metadata = json.load(open(prefix + '.json'))
    data = np.load(prefix + '.npz')
    Rc, R90 = data['Rc'], data['R90']
    for x, y in zip(Rc, R90):
        ax.plot(x, y, '-', lw=2, c=color, label="_nolabel_", alpha=0.3)
    ax.plot(Rc[0], R90[0], '-', lw=3, c=color,
            label=metadata['shape_label'], alpha=0.8)

ax.legend(title=metadata['wave_label'])

ax.set(
    yscale='linear',
    xscale='linear',
    xlim=[0.0, 5.1],
    ylim=[0.0, 5.1],
    xlabel=r"Projected dimensionless radius of curvature: $\widetilde{R}_{c}{}'$",
    ylabel=r"Projected dimensionless perpendicular radius: $\widetilde{R}_{90}{}'$",
)        

sns.despine()
fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
