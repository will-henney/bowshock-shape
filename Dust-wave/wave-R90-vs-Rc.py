import sys
import json
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns

try:
    WAVE_ID = sys.argv[1]
    BASE_SHAPE_IDS = sys.argv[2:]
except:
    sys.exit(f"Usage: {sys.argv[0]} WAVE_ID SHAPE_ID [SHAPE_ID ...]")

plotfile = sys.argv[0].replace('.py', f'-{WAVE_ID}.pdf')

range_ = [0.35, 9.5]
ticks = [0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0]

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
ax.plot([5./3.], [np.sqrt(3.0)], '+', c='w', ms=10, alpha=1.0, zorder=100)

nshapes = len(BASE_SHAPE_IDS)
colors = sns.color_palette(n_colors=nshapes)
for base_shape_id, color in zip(BASE_SHAPE_IDS, colors):
    prefix = f"{base_shape_id}-wave-{WAVE_ID}"
    metadata = json.load(open(prefix + '.json'))
    data = np.load(prefix + '.npz')
    Rc, R90 = data['Rc'], data['R90']
    for x, y in zip(Rc, R90):
        ax.plot(x, y, '.', ms=5, mec="none", c=color, label="_nolabel_", alpha=0.1)
    ax.plot([], [], '.', ms=9, mec="none", c=color,
            label=metadata['shape_label'], alpha=0.8)

ax.legend(title=metadata['wave_label'])

ax.set(
    xlim=range_, ylim=range_,
    xlabel=r'Planitude: $\Pi = R_c/R_0$',
    ylabel=r'Alatude: $\Lambda = R_{90}/R_0$',
    xscale='log', yscale='log')
ax.xaxis.set_major_locator(mticker.FixedLocator(ticks))
ax.yaxis.set_major_locator(mticker.FixedLocator(ticks))
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%0g'))
ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%0g'))

sns.despine()
fig.tight_layout(pad=1.0)
fig.savefig(plotfile)
print(plotfile, end='')
