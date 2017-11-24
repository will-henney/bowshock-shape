T=[["Source", "R0/100", "Rc/R0", "R90/R0", "Rm90/R0"], ["A16-12mic-i00", "0.44", "2.59", "1.69", "1.69"], ["A16-12mic-i30", "0.42", "2.90", "1.97", "1.88"], ["A16-12mic-i60", "0.54", "1.90", "1.45", "1.42"], ["A16-20mic-i00", "0.41", "2.05", "1.75", "1.87"], ["A16-20mic-i30", "0.40", "2.41", "1.96", "1.97"], ["A16-20mic-i60", "0.43", "4.94", "1.78", "1.82"], ["A16-Halpha-i00", "0.52", "1.83", "1.57", "1.79"], ["A16-Halpha-i30", "0.63", "2.25", "1.53", "1.54"], ["A16-Halpha-i60", "0.81", "1.37", "1.38", "1.38"], ["R97-Halpha-i00", "0.92", "1.55", "2.16", "2.09"]]
import sys
import numpy as np
from scipy.interpolate import interp1d
import statsmodels.api as sm
from astropy.table import Table
from matplotlib import pyplot as plt
import matplotlib.ticker
import seaborn as sns
sys.path.append("../Dust-wave")
import bow_projection as bp
import bow_diagnostic


plotfile = sys.argv[0].replace('.py', '.pdf')

table = Table(rows=T[1:], names=T[0], dtype=[str] + [float]*4)

# Take average +/- std of the +ve and -ve R90
R90stack = np.stack([table['R90/R0'], table['Rm90/R0']])
table['R90'] = np.nanmean(R90stack, axis=0)
table['dR90'] = np.nanstd(R90stack, axis=0)
table.remove_columns(['R90/R0', 'Rm90/R0'])

sns.set_style('ticks')
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
ax.plot([0.0, 10.0], [0.0, 10.0], lw=0.5, alpha=0.5, color='k', zorder=-1)

# Put a cross at the Wilkinoid coordinates: [5/3, sqrt(3)]
ax.plot([5./3.], [np.sqrt(3.0)], '+', c='w', ms=10, alpha=1.0)
# And plot the projected wilkinoids 
bp.N_NEIGHBORHOOD = 50
bp.DEGREE_POLY_NEIGHBORHOOD = 2
bp.SCALE_NEIGHBORHOOD = 0.03
bp.DEGREE_POLY_NEIGHBORHOOD_90 = 2
bp.SCALE_NEIGHBORHOOD_90 = 0.01
shape = bp.wilkinoid_R_theta
th_inf = bp.theta_infinity(shape)
inc = np.linspace(0.0, th_inf - np.pi/2, 50)
tab = bow_diagnostic.parameter_table(inc, shape)
Rc, R90 = tab['tilde R_c prime'], tab['tilde R_90 prime']
ax.plot(Rc, R90, '-', c='w', label="_nolabel_", lw=0.6, alpha=0.9)
sini = (0.5 + np.arange(20))/20
inc_e = np.arcsin(sini)
tab_e = bow_diagnostic.parameter_table(inc_e, shape)
Rc_e, R90_e = tab_e['tilde R_c prime'], tab_e['tilde R_90 prime']
ax.scatter(Rc_e, R90_e, marker='|', s=3**2,
           linewidths=0.1, edgecolors='none',
           c='w', alpha=0.5, label="_nolabel_")


models = ["A16-12mic", "A16-20mic", "A16-Halpha", "R97-Halpha"]

colors = sns.color_palette(n_colors=len(models))
for model, color in zip(models, colors):
    mask = [s.startswith(model) for s in table['Source']]
    data = table[mask]
    ax.plot(data['Rc/R0'], data['R90'], '-', c=color, label=model, lw=1.5, alpha=0.9)
    # Put a dot at the i=0 case
    ax.plot(data['Rc/R0'][0:1], data['R90'][0:1], 'o', mec='none', c=color, label="_nolabel_", alpha=0.7)



ax.legend(ncol=1, fontsize='small', title='Simulations',
          frameon=True, loc="lower right")
ax.set(
    xlim=[0.0, 5.1],
    ylim=[0.0, 5.1],
    yticks=range(6),
#    ylim=[-3.0, 1.1],
    xlabel=r"Projected planitude: $\Pi'$",
    ylabel=r"Projected alatude: $\Lambda'$",
)        

sns.despine()
fig.tight_layout(pad=0.5)
fig.savefig(plotfile)
print(plotfile, end='')
