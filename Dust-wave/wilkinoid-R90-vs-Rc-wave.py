import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker
import seaborn as sns
import bow_projection as bp
import bow_diagnostic
import dragoid_shape
import ancantoid_shape
import standing_wave

bp.N_NEIGHBORHOOD = 50
bp.DEGREE_POLY_NEIGHBORHOOD = 2
bp.SCALE_NEIGHBORHOOD = 0.03
bp.DEGREE_POLY_NEIGHBORHOOD_90 = 2
bp.SCALE_NEIGHBORHOOD_90 = 0.01

try:
    amplitude = float(sys.argv[1])
    wavenumber = float(sys.argv[2])
except:
    sys.exit(f"Usage: {sys.argv[0]} AMPLITUDE WAVENUMBER")

plotfile = sys.argv[0].replace(
    '.py', f'-A{int(100*amplitude):03d}-N{int(10*wavenumber):02d}.pdf')

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
# And plot the projected wilkinoids 
shape = bp.wilkinoid_R_theta
th_inf = bp.theta_infinity(shape)
inc = np.linspace(0.0, th_inf - np.pi/2, 50)
tab = bow_diagnostic.parameter_table(inc, shape)
Rc, R90 = tab['tilde R_c prime'], tab['tilde R_90 prime']
ax.plot(Rc, R90, '-', c='w', label="_nolabel_", lw=0.6, alpha=0.9)
ax.plot(Rc, R90, 'x', ms=0.2, c='w', label="_nolabel_", alpha=0.5)

base_shapes = [
    bp.paraboloid_R_theta,
    bp.wilkinoid_R_theta, 
    dragoid_shape.Dragoid(alpha=1.0),
    bp.Spline_R_theta_from_function(ngrid=1000,
                                    shape_func=bp.cantoid_R_theta,
                                    shape_func_pars=(0.1,)),
    ancantoid_shape.Ancantoid(xi=0.8, beta=0.005, n=301),
]

shape_labels = [
    "Paraboloid",
    "Wilkinoid",
    r"Dragoid $\alpha_\mathrm{drag} = 1.0$",
    r"Cantoid $\beta = 0.1$",
    r"Ancantoid $\xi = 0.8$, $\beta = 0.05$" 
]

nshapes = len(base_shapes)
cols = sns.color_palette(n_colors=nshapes)

shapes = [standing_wave.StandingWave(base_shape,
                                     amplitude=amplitude,
                                     wavenumber=wavenumber)
          for base_shape in base_shapes]


nphases = 21
ninc = 50
# nphases = 5
# ninc = 20
phases = np.linspace(0.0, 0.5, nphases)
mus = np.linspace(0.0, 1.0, ninc)
inclinations = np.arcsin(mus)
for phase in phases[::-1]:
    for shape, shape_label, col in zip(shapes, shape_labels, cols):
        shape.phase = phase
        tab = bow_diagnostic.parameter_table(inclinations, shape)
        Rc, R90 = tab['tilde R_c prime'], tab['tilde R_90 prime']
        ax.plot(Rc, R90, '-', lw=2, c=col,
                label="_nolabel_", alpha=0.3)
        # And plot the inc=0 case bolder
        label = shape_label if phase == 0.0 else "_nolabel_"
        ax.plot(Rc[0:1], R90[0:1], '-', lw=3, c=col,
                label=label, alpha=0.8)

ax.legend(title= f"Perturbation $A = {amplitude:.03f}$, $N = {wavenumber:0.1f}$")

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
