import sys
import numpy as np
from scipy.interpolate import interp1d
import statsmodels.api as sm
from astropy.table import Table
from matplotlib import pyplot as plt
import matplotlib.ticker
import seaborn as sns
from simulation_shape import Simulation
sys.path.append("../Dust-wave")
import bow_projection as bp
import bow_diagnostic

plotfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('ticks')
fig, axes = plt.subplots(3, 1, figsize=(4, 5))


bp.N_NEIGHBORHOOD = 200
bp.DEGREE_POLY_NEIGHBORHOOD = 1
bp.SCALE_NEIGHBORHOOD = (60./180.) # => +/-60 deg at i=0
bp.DEGREE_POLY_NEIGHBORHOOD_90 = 2
bp.SCALE_NEIGHBORHOOD_90 = 0.01

models = ["M17-HD2040", "M17-MHD2040-AllB7"][::-1]
labels = ["HD", "MHD"][::-1]

colors = sns.color_palette(n_colors=len(models))[::-1]

sini = np.linspace(0.0, 1.0, 100)
inc = np.arcsin(sini)
R0p, Rc, R90 = [], [], []
for model in models:
    # First do the "theoretical" tracks
    shape = Simulation(name=model, force_open=True, cheby_degree=12)
    tab = bow_diagnostic.parameter_table(inc, shape)
    R0p.append(tab['R_0 prime'])
    Rc.append(tab['tilde R_c prime'])
    R90.append(tab['tilde R_90 prime'])

axes[0].hist(Rc,  density=True, bins=12, alpha=0.7, range=[0.0, 6.0], color=colors, label=labels)
axes[1].hist(R90, density=True, bins=12, alpha=0.7, range=[1.0, 3.0], color=colors, label=labels)
axes[2].hist(R0p, density=True, bins=12, alpha=0.7, range=[0.0, 3.0], color=colors, label=labels)

axes[1].legend(ncol=1, fontsize='small',
             title='Simulation\ncontact discontinuity\n(Meyer et al. 2017)',
             frameon=True).get_title().set_size('small')
axes[0].set(
    xlabel=r"Projected planitude: $\Pi'$",
    yticks=[0, 1],
)
axes[1].set(
    xlabel=r"Projected alatude: $\Lambda'$",
    ylabel="Probability density",
)
axes[2].set(
    xlabel=r"Projected apex distance: $R_{0}'/R_{0}$",
)

sns.despine(trim=True)
fig.tight_layout(pad=0.5)
fig.savefig(plotfile)
print(plotfile, end='')
