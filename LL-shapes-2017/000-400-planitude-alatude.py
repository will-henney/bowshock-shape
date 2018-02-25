import sys
import numpy as np
from scipy.interpolate import interp1d
import statsmodels.api as sm
from astropy.table import Table
from matplotlib import pyplot as plt
import matplotlib.ticker
import seaborn as sns


plotfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('ticks')
sns.set_color_codes()
fig, ax = plt.subplots(figsize=(4, 4))

Rc_grid = np.linspace(0.0, 10.0, 2000)
R90_T0_grid = np.sqrt(2*Rc_grid)
R90_T1_grid = np.sqrt(2*Rc_grid - 1.0)
R90_T1_grid[~np.isfinite(R90_T1_grid)] = 0.0 

ax.fill_between(Rc_grid, R90_T1_grid, R90_T0_grid, color='k', alpha=0.2)
ax.fill_between(Rc_grid, R90_T0_grid, color='k', alpha=0.1)
ax.plot(Rc_grid, R90_T0_grid, c='k', lw=0.5)
ax.axhline(1.0, lw=0.5, alpha=0.5, color='k', zorder=-1)
ax.axvline(1.0, lw=0.5, alpha=0.5, color='k', zorder=-1)
ax.axhline(2.0, ls=':', lw=0.5, alpha=0.5, color='k', zorder=-1)
ax.axvline(2.0, ls=':', lw=0.5, alpha=0.5, color='k', zorder=-1)
ax.plot([0.9, 10.0], [0.9, 10.0], lw=0.5, alpha=0.5, color='k', zorder=-1)

# Put a cross at the Wilkinoid coordinates: [5/3, sqrt(3)]
ax.plot([5./3.], [np.sqrt(3.0)], '+', c='w', ms=10, alpha=1.0)

labels = {
    "Inner": 1.7,
    "Outer": 2.7,
    "Ridge A": 3.4,
}
for arcid, Pi, Lam_minus, Lam_plus in [
        ["Inner", 2.692, 2.795, 2.021],
        ["Ridge A", 3.664, 3.225, 2.242],
        ["Ridge B", 3.318, 3.162, 2.257],
        ["Outer", 2.186, 2.368, 2.512],
        ["Inner 45", 2.702, 2.518, 2.111],
        ["Ridge 45", 4.441, 2.950, 2.846],
        ["Outer 45", 3.111, 2.460, 2.582],
]:

    if "Ridge" in arcid:
        color = 'r'
    elif "Inner" in arcid:
        color = 'm'
    elif "Outer" in arcid:
        color = 'c'
    if "45" in arcid:
        sym = '.'
        lw = 0.7
    else:
        sym = 'o'
        lw = 1.5
    ax.plot([Pi, Pi], [Lam_minus, Lam_plus], '-', color=color, lw=lw)
    ax.plot([Pi, Pi], [Lam_minus, Lam_plus], sym, color=color)
    if arcid in labels:
        ax.text(Pi, labels[arcid], arcid.split()[0], color=color, ha='center')

ax.set(
    xlim=[0., 5.1],
    ylim=[0., 5.1],
    #yticks=range(6),
    xlabel=r"Projected planitude: $\Pi'$",
    ylabel=r"Projected alatude: $\Lambda'$",
)        
sns.despine()
fig.tight_layout(pad=0.5)
fig.savefig(plotfile)
print(plotfile, end='')
