import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker
import seaborn as sns
import bow_projection as bp
import ancantoid_shape
import bow_diagnostic

try: 
    xiset = sys.argv[1]
    plotfile = sys.argv[0].replace('.py', f'-{xiset}.pdf')
    assert xiset in 'ab'
    istart = -2 if xiset == 'a' else -1
except:
    sys.exit(f"Usage: {sys.argv[0]} a|b")

sns.set_style('ticks')
fig, ax = plt.subplots(figsize=(5, 5))

bp.N_NEIGHBORHOOD = 50
bp.DEGREE_POLY_NEIGHBORHOOD = 2
bp.SCALE_NEIGHBORHOOD = 0.03
bp.DEGREE_POLY_NEIGHBORHOOD_90 = 2
bp.SCALE_NEIGHBORHOOD_90 = 0.01

left_annotate_pars = dict(xytext=(-5, 5), ha='right', va='bottom')
right_annotate_pars = dict(xytext=(5, -5), ha='left', va='top')


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

XI_LIST = [None, 1.0, 0.8, 0.4]
BETA_LIST = [0.2, 0.1, 0.05, 0.02, 0.005, 1e-3]
nxi, nbeta = len(XI_LIST), len(BETA_LIST)

# Put a cross at the Wilkinoid coordinates: [5/3, sqrt(3)]
ax.plot([5./3.], [np.sqrt(3.0)], '+', c='w', ms=10, alpha=1.0)
# And plot the projected wilkinoids 
shape = bp.wilkinoid_R_theta
th_inf = bp.theta_infinity(shape)
inc = np.linspace(0.0, th_inf - np.pi/2, 50)
tab = bow_diagnostic.parameter_table(inc, shape)
Rc, R90 = tab['tilde R_c prime'], tab['tilde R_90 prime']
ax.plot(Rc, R90, '-', c='w', label="_nolabel_", lw=0.6, alpha=0.9)
ax.plot(Rc, R90, 'x', ms=0.2, c='w', label="_nolabel_", alpha=0.5)


cols = sns.color_palette('magma', n_colors=nxi)
annot_pars_list = [right_annotate_pars]*2 + [left_annotate_pars]*2 
for beta in BETA_LIST[::-1]:
    for xi, col, annot_pars in list(zip(XI_LIST, cols, annot_pars_list))[istart::-2]:

        if beta == BETA_LIST[0]:
            label = "Cantoid" if xi is None else fr"Ancantoid $\xi = {xi:.1f}$"
        else:
            label = "_nolabel_"

        if xi is None:
            shape = bp.Spline_R_theta_from_function(
                ngrid=1000,
                shape_func=bp.cantoid_R_theta,
                shape_func_pars=(beta,))
        else:
            shape = ancantoid_shape.Ancantoid(xi=xi, beta=beta, n=301)

        th_inf = bp.theta_infinity(shape)
        inc = np.linspace(0.0, th_inf - np.pi/2, 50)
        tab = bow_diagnostic.parameter_table(inc, shape)
        Rc, R90 = tab['tilde R_c prime'], tab['tilde R_90 prime']

        ax.plot(Rc, R90, '-', c=col, label=label, lw=0.4, alpha=0.9)
        ax.plot(Rc, R90, 'x', ms=0.2, c=col, label="_nolabel_", alpha=0.5)
        # Put a dot at the i=0 case
        ax.plot(Rc[0:1], R90[0:1], '.', mec='none', c=col, label="_nolabel_", alpha=0.7)
        # Label the dot with the cross-over inclination
        ax.annotate(rf'$\beta = \mathrm{{{beta:g}}}$',
                    xy=(Rc[0], R90[0]),
                    textcoords='offset points',
                    fontsize='x-small', color=col, **annot_pars)


ax.legend(ncol=1, fontsize='small', frameon=True)
ax.set(
    yscale='linear',
    xscale='linear',
    xlim=[0.0, 5.1],
    ylim=[0.0, 5.1],
#    ylim=[-3.0, 1.1],
    xlabel=r"Projected dimensionless radius of curvature: $\widetilde{R}_{c}{}'$",
    ylabel=r"Projected dimensionless perpendicular radius: $\widetilde{R}_{90}{}'$",
)        

sns.despine()
fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
