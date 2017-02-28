T=[["073-227", "147.26", "0.96", "6.47", "-999.00", "3.12", "1.61", "3.20", "-999.00", "5.22"], ["109-246", "89.67", "1.21", "27.23", "1.83", "-999.00", "1.90", "6.22", "2.65", "-999.00"], ["000-400", "254.03", "0.79", "3.37", "2.46", "2.34", "1.46", "2.21", "2.59", "2.25"], ["005-514", "262.64", "1.14", "2.76", "-999.00", "-999.00", "1.61", "1.71", "1.50", "-999.00"], ["012-407", "231.47", "-999.00", "-999.00", "-999.00", "-999.00", "2.29", "1.86", "-999.00", "1.92"], ["030-524", "234.09", "0.37", "1.99", "2.33", "2.91", "0.63", "3.84", "2.17", "3.10"], ["042-628", "259.60", "1.78", "2.77", "1.67", "-999.00", "3.08", "1.87", "1.58", "2.32"], ["LL1", "198.63", "1.90", "2.61", "3.05", "2.05", "3.06", "2.46", "2.14", "2.28"], ["069-601", "212.20", "0.45", "2.90", "2.42", "3.39", "0.89", "2.05", "2.13", "2.09"], ["4285-458", "721.18", "-999.00", "-999.00", "-999.00", "-999.00", "1.91", "2.56", "1.68", "2.44"], ["LL3", "566.33", "1.50", "2.35", "-999.00", "-999.00", "3.12", "1.94", "1.54", "2.41"], ["LL4", "593.14", "1.42", "5.81", "4.67", "3.35", "2.54", "2.35", "2.69", "2.93"], ["4468-605", "471.31", "1.31", "2.47", "1.35", "1.52", "2.50", "1.80", "2.58", "1.58"], ["116-3101", "463.92", "1.00", "1.41", "1.72", "1.74", "1.46", "1.42", "1.52", "1.84"], ["266-558", "218.16", "1.12", "1.79", "-999.00", "1.25", "2.05", "1.92", "1.98", "2.94"], ["308-3036", "484.14", "1.40", "1.35", "-999.00", "1.05", "2.51", "1.50", "1.49", "2.08"], ["LL5", "369.55", "1.46", "2.63", "1.59", "2.69", "2.89", "3.03", "2.22", "4.40"], ["LL6", "485.83", "1.63", "6.01", "3.77", "4.31", "3.72", "5.57", "-999.00", "2.75"]]
import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker
import seaborn as sns
from astropy.table import Table

cols = ['Source', 'D',
        'R0 in', 'Rc in', 'R90 in', 'Rm90 in',
        'R0 out', 'Rc out', 'R90 out', 'Rm90 out',
]
dtype = [str] + [float]*9
table = Table(rows=T, names=cols, dtype=dtype)

Rcols = [_ for _ in cols if _.startswith('R')]
for col in Rcols:
    m = table[col] < 0.0
    table[col][m] = np.nan

# Take average +/- std of the +ve and -ve R90
R90stack = np.stack([table['R90 out'], table['Rm90 out']])
table['R90'] = np.nanmean(R90stack, axis=0)
table['dR90'] = np.nanstd(R90stack, axis=0)
table.remove_columns(['R90 out', 'Rm90 out'])

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

ax.scatter(table['Rc out'], table['R90'])
ax.errorbar(table['Rc out'], table['R90'], yerr=table['dR90'], fmt='none', alpha=0.3)

ax.set(
    yscale='linear',
    xscale='linear',
    xlim=[0.0, 8.1],
    ylim=[0.0, 8.1],
#    ylim=[-3.0, 1.1],
    xlabel=r"Projected dimensionless radius of curvature: $\widetilde{R}_{c}{}'$",
    ylabel=r"Projected dimensionless perpendicular radius: $\widetilde{R}_{90}{}'$",
)        


fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
