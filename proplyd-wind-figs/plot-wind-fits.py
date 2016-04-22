"""Plot data from table of pressures and fluxes from proplyd arc fits
"""
from astropy.table import Table
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns

AU = 1.49597870691e13
PC = 3.085677582e18
d_Orion = 440.0
k_Boltzmann = 1.3806503e-16
cos80 = 0.173648177667

tab = Table.read('../doc/wind-fits.tab', format='ascii.tab')

sources = sorted(set(tab['Source']))
n = len(sources)
colors = sns.color_palette('Set1', n)

# Physical separations in parsec
D = tab['D'] * d_Orion*AU/PC
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.loglog(D, tab['F(star)'], c='k', alpha=0.1, lw=10, label='')
ax2.loglog(D, tab['P(wind)']/k_Boltzmann, c='k', alpha=0.1, lw=10, label='')
mm = (tab['F(ph)/F(*)'] > 0.707) & (tab['F(ph)/F(*)'] < 1.414)
for source, color in zip(sources, colors):
    m = tab['Source'] == source
    Dprime = tab['D\''][m][0] * d_Orion*AU/PC
    F = tab['F(photo)'][m][0]
    ax1.loglog([Dprime, Dprime/cos80], [F, F], ':', c=color, alpha=0.4, label='')
    ax1.loglog(D[m], tab['F(photo)'][m],
               '-', c=color, alpha=0.4, label='')
    ax1.loglog(D[m & mm], tab['F(photo)'][m & mm],
               'o-', lw=3, c=color, label=source)
    ax2.loglog(D[m], tab['P(in)'][m]/k_Boltzmann,
               '-', c=color, alpha=0.4, label='')
    ax2.loglog(D[m & mm], tab['P(in)'][m & mm]/k_Boltzmann,
               'o-', c=color, lw=3, label=source)
ax2.legend(ncol=2, loc='lower left')
ax2.set_xlim(0.008, 0.3)
ax2.set_ylim(1e7, 4e9)
ax1.set_ylim(2e12, 3e14)
ax1.set_ylabel(r'Ionizing Flux, $\mathrm{cm^{-2}\ s^{-1}}$')
ax2.set_ylabel(r'Stagnation Pressure: $P/k$, $\mathrm{cm^{-3}\ K}$')
ax2.set_xlabel('Distance, parsec')
fig.set_size_inches(5, 8)
fig.tight_layout()
fig.savefig('plot-wind-fits.pdf')


