import sys
import numpy as np
from scipy.interpolate import interp1d
import statsmodels.api as sm
from matplotlib import pyplot as plt
import matplotlib.ticker
import seaborn as sns
import bow_projection as bp
import dragoid_shape
import bow_diagnostic


plotfile = sys.argv[0].replace('.py', '.pdf')

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
sini = np.linspace(0.0, 1.0, 20)
inc_e = np.arcsin(sini)
tab_e = bow_diagnostic.parameter_table(inc_e, shape)
Rc_e, R90_e = tab_e['tilde R_c prime'], tab_e['tilde R_90 prime']
ax.scatter(Rc_e, R90_e, marker='|', s=3**2,
           linewidths=0.1, edgecolors='none',
           c='w', alpha=0.5, label="_nolabel_")

# For dragoids, widen the region for fitting R_c
bp.SCALE_NEIGHBORHOOD = 0.08

ALPHA_LIST = [0.25, 1.0, 1.0, 1.0, 4.0, 4.0]
MU_LIST = [None, None, 0.05, 0.2, 0.2, 0.8]
nalpha = len(ALPHA_LIST)
cols = sns.color_palette(n_colors=nalpha)
for alpha, mu, col in list(zip(ALPHA_LIST, MU_LIST, cols))[::-1]:
    shape = dragoid_shape.Dragoid(alpha=alpha, mu=mu)
    th_inf = bp.theta_infinity(shape)
    inc = np.linspace(0.0, th_inf - np.pi/2, 50)
    tab = bow_diagnostic.parameter_table(inc, shape)
    Rc, R90 = tab['tilde R_c prime'], tab['tilde R_90 prime']

    Rcs = sm.nonparametric.lowess(Rc, inc, frac=0.3,
                                  is_sorted=True, return_sorted=False)
    # Rcs = Rc
    Rcs[0] = Rc[0]
    ax.plot(Rcs, R90, '-', c=col, label=shape.label, lw=0.7, alpha=0.9)

    # Get points evenly spaced in sin i
    sini = np.linspace(0.0, 1.0, 20)
    inc_e = np.arcsin(sini)
    inc_e = inc_e[inc_e < th_inf - np.pi/2]
    # Interpolate to get the even probability points
    Rc_e = interp1d(inc, Rcs)(inc_e)
    R90_e = interp1d(inc, R90)(inc_e)

    ax.scatter(Rc_e, R90_e, marker='|', s=3**2,
               linewidths=0.1, edgecolors='none',
               c=col, alpha=0.5, label="_nolabel_")

    # Put a dot at the i=0 case
    ax.plot(Rc[0:1], R90[0:1], 'o', mec='none', c=col, label="_nolabel_", alpha=0.7)


ax.legend(ncol=1, fontsize='small', title='Dragoids',
          frameon=True, loc="lower right")
ax.set(
    yscale='linear',
    xscale='linear',
    xlim=[0.9, 2.7],
    ylim=[0.9, 2.7],
#    ylim=[-3.0, 1.1],
    xlabel=r"Projected planitude: $\Pi'$",
    ylabel=r"Projected alatude: $\Lambda'$",
)        

sns.despine()
fig.tight_layout(pad=0.5)
fig.savefig(plotfile)
print(plotfile, end='')
