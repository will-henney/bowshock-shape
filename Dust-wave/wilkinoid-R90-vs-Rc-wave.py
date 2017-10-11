import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker
import seaborn as sns
import bow_projection as bp
import bow_diagnostic
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


nphases = 11
ninc = 50
phases = np.linspace(0.0, 0.5, nphases)
mus = np.linspace(0.0, 1.0, ninc)
inclinations = np.arcsin(mus)
shape = standing_wave.StandingWave(
    bp.wilkinoid_R_theta,
    amplitude=amplitude, wavenumber=wavenumber)
cols = sns.color_palette('magma', n_colors=nphases)
for phase, col in list(zip(phases, cols))[::-1]:
    shape.phase = phase
    th_inf = bp.theta_infinity(shape)
    tab = bow_diagnostic.parameter_table(inclinations, shape)
    Rc, R90 = tab['tilde R_c prime'], tab['tilde R_90 prime']
    ax.plot(Rc, R90, '.',  mec='none', c=col, label="_nolabel_", alpha=0.8)

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
