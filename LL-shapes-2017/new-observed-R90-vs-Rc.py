T=[["073-227", "147.26", "0.65", "6.47", "-999.00", "3.12", "1.09", "3.20", "-999.00", "5.22"], ["109-246", "89.67", "1.34", "27.23", "1.83", "-999.00", "2.11", "6.22", "2.65", "-999.00"], ["000-400", "254.03", "0.31", "3.37", "2.46", "2.34", "0.58", "2.21", "2.59", "2.25"], ["005-514", "262.64", "0.44", "2.76", "-999.00", "-999.00", "0.61", "1.71", "1.50", "-999.00"], ["012-407", "231.47", "-999.00", "-999.00", "-999.00", "-999.00", "0.99", "1.86", "-999.00", "1.92"], ["030-524", "234.09", "0.16", "1.99", "2.33", "2.91", "0.27", "3.84", "2.17", "3.10"], ["042-628", "259.60", "0.69", "2.77", "1.67", "-999.00", "1.19", "1.87", "1.58", "2.32"], ["LL1", "198.63", "0.96", "2.61", "3.05", "2.05", "1.54", "2.46", "2.14", "2.28"], ["069-601", "212.20", "0.21", "2.90", "2.42", "3.39", "0.42", "2.05", "2.13", "2.09"], ["4285-458", "721.18", "-999.00", "-999.00", "-999.00", "-999.00", "0.27", "2.56", "1.68", "2.44"], ["LL3", "566.33", "0.26", "2.35", "-999.00", "-999.00", "0.55", "1.94", "1.54", "2.41"], ["LL4", "593.14", "0.24", "5.81", "4.67", "3.35", "0.43", "2.35", "2.69", "2.93"], ["4468-605", "471.31", "0.28", "2.47", "1.35", "1.52", "0.53", "1.80", "2.58", "1.58"], ["116-3101", "463.92", "0.22", "1.41", "1.72", "1.74", "0.31", "1.42", "1.52", "1.84"], ["266-558", "218.16", "0.51", "1.79", "-999.00", "1.25", "0.94", "1.92", "1.98", "2.94"], ["308-3036", "484.14", "0.29", "1.35", "-999.00", "1.05", "0.52", "1.50", "1.49", "2.08"], ["LL5", "369.55", "0.40", "2.63", "1.59", "2.69", "0.78", "3.03", "2.22", "4.40"], ["LL6", "485.83", "0.33", "6.01", "3.77", "4.31", "0.77", "5.57", "-999.00", "2.75"]]
TRSG=[["alphaori", "3.70", "1.43", "1.38", "1.43"], ["uuaur", "0.78", "1.31", "1.22", "1.39"], ["rleo", "0.85", "1.43", "1.29", "1.39"], ["rhya", "0.89", "1.48", "1.46", "1.69"], ["v1943sgr", "0.53", "1.34", "1.42", "1.20"], ["xpav", "0.48", "1.51", "1.48", "1.50"], ["mucep", "0.70", "1.56", "1.52", "1.26"], ["cwleo", "4.68", "1.53", "1.41", "1.60"], ["epaqr", "0.31", "2.18", "1.62", "-999.00"], ["khicyg", "2.86", "1.47", "1.35", "-999.00"], ["rcas", "0.85", "1.65", "1.56", "1.46"], ["rtvir", "0.87", "1.13", "1.25", "1.38"], ["waql", "0.42", "1.28", "1.39", "-999.00"], ["wpic", "0.28", "1.60", "1.45", "-999.00"], ["rscl", "0.48", "1.23", "1.31", "-999.00"], ["tetaps", "0.57", "1.53", "-999.00", "1.29"]]
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
cols2 = ['Source',
        'R0 out', 'Rc out', 'R90 out', 'Rm90 out',
]
dtype = [str] + [float]*9
dtype2 = [str] + [float]*4
table_LL = Table(rows=T, names=cols, dtype=dtype)
table_RSG = Table(rows=TRSG, names=cols2, dtype=dtype2)

Rcols = [_ for _ in cols if _.startswith('R')]

for table in table_LL, table_RSG:
    for col in Rcols:
        if col in table.colnames:
            m = table[col] < 0.0
            table[col][m] = np.nan

    # Take average +/- std of the +ve and -ve R90
    R90stack = np.stack([table['R90 out'], table['Rm90 out']])
    table['R90'] = np.nanmean(R90stack, axis=0)
    table['dR90'] = np.nanstd(R90stack, axis=0)
    table.remove_columns(['R90 out', 'Rm90 out'])

table_LL.write('ll-arcs-radii.tab', overwrite=True, format='ascii.tab')
table_RSG.write('rsg-arcs-radii.tab', overwrite=True, format='ascii.tab')

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

ax.scatter(table_LL['Rc out'], table_LL['R90'], s=40*table_LL['R0 out'])
ax.errorbar(table_LL['Rc out'], table_LL['R90'], yerr=table_LL['dR90'], fmt='none', alpha=0.3)

ax.scatter(table_RSG['Rc out'], table_RSG['R90'], s=10, c='r', alpha=0.8)
ax.errorbar(table_RSG['Rc out'], table_RSG['R90'], yerr=table_RSG['dR90'], fmt='none', alpha=0.3)

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
