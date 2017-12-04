T=[["Source", "R0/pc", "Rc/R0", "R90/R0", "Rm90/R0"], ["M17-MHD2040-AllB7", "0.38", "3.57", "1.95", "1.92"], ["M17-MHD2040-AllB7-Halpha-i00", "0.40", "2.83", "1.90", "1.89"], ["M17-MHD2040-AllB7-60mic-i30", "0.48", "2.14", "1.74", "1.73"], ["M17-MHD2040-AllB7-60mic-i45", "0.56", "1.70", "1.59", "1.61"], ["M17-MHD2040-AllB7-60mic-i60", "0.67", "1.51", "1.47", "1.46"], ["M17-HD2040", "0.66", "1.78", "1.75", "1.75"], ["M17-HD2040-Halpha-i00", "0.71", "2.33", "1.83", "1.85"], ["M17-HD2040-Halpha-i00-CD", "0.59", "1.68", "1.93", "1.93"], ["M17-HD2040-Halpha-i00-BS", "0.75", "1.91", "1.83", "1.83"], ["M17-HD2040-60mic-i30", "1.17", "1.89", "1.80", "1.79"], ["M17-HD2040-60mic-i45", "1.36", "1.87", "1.90", "1.59"], ["M17-HD2040-60mic-i60", "1.70", "1.66", "1.55", "1.59"], ["M17-MHD2040-AllB7-60mic-i90", 1.05, -999.0, -999.0, -999.0], ["M17-HD2040-60mic-i90", 3.0, -999.0, -999.0, -999.0]]
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
fig, ax = plt.subplots(figsize=(4, 4))


table = Table(rows=T[1:], names=T[0], dtype=[str] + [float]*4)

# Take average +/- std of the +ve and -ve R90
R90stack = np.stack([table['R90/R0'], table['Rm90/R0']])
table['R90'] = np.nanmean(R90stack, axis=0)
table['dR90'] = np.nanstd(R90stack, axis=0)
table.remove_columns(['R90/R0', 'Rm90/R0'])


def select_marker_style(s):
    if 'Halpha' in s:
        return '^'
    elif '60mic' in s:
        return 's'
    else:
        return 'o'

def select_marker_size(s):
    if 'i00' in s:
        return 5
    elif 'i30' in s:
        return 6
    elif 'i45' in s:
        return 5
    elif 'i60' in s:
        return 4
    elif 'i90' in s:
        return 3
    else:
        return 5

def select_inclination(s):
    if 'i00' in s:
        return 0.0
    elif 'i30' in s:
        return 30.0
    elif 'i45' in s:
        return 45.0
    elif 'i60' in s:
        return 60.0
    elif 'i90' in s:
        return 90.0
    else:
        return 0.0


table['marker style'] = [select_marker_style(s) for s in table['Source']]
table['marker size'] = [select_marker_size(s) for s in table['Source']]
table['inclination'] = [select_inclination(s) for s in table['Source']]

bp.N_NEIGHBORHOOD = 200
bp.DEGREE_POLY_NEIGHBORHOOD = 1
bp.SCALE_NEIGHBORHOOD = (60./180.) # => +/-60 deg at i=0
bp.DEGREE_POLY_NEIGHBORHOOD_90 = 2
bp.SCALE_NEIGHBORHOOD_90 = 0.01

models = ["M17-HD2040", "M17-MHD2040-AllB7"][::-1]
labels = ["HD", "MHD"][::-1]

colors = sns.color_palette(n_colors=len(models))[::-1]
for model, label, color in zip(models, labels, colors):
    R0 = 0.378 if 'MHD' in model else 0.661
    ax.axhline(R0, ls=':', lw=0.5, color=color)

    inc = np.linspace(0.0, np.pi/2, 100)

    # First do the "theoretical" tracks
    shape = Simulation(name=model, force_open=True, cheby_degree=12)
    tab = bow_diagnostic.parameter_table(inc, shape)
    R0p = tab['R_0 prime']
    ax.plot(np.degrees(inc), R0*R0p, c=color, label=label, lw=2.5, alpha=0.5)

    # Repeat for the closed shape, but with a thin line
    shape = Simulation(name=model, force_open=False, cheby_degree=12)
    tab = bow_diagnostic.parameter_table(inc, shape)
    R0pp = tab['R_0 prime']
    ax.plot(np.degrees(inc), R0*R0pp, '-', c=color, label="_nolabel_", lw=0.5, alpha=1.0)

    mask = [s.startswith(model) and not (s.endswith('BS') or s.endswith('CD')) 
            for s in table['Source']]
    data = table[mask]

    for row in data:
        FIX = 0.63 if label == "HD" and "60mic" in row["Source"] else 1.0
        ax.scatter(row['inclination'], row['R0/pc']*FIX,
                   marker=row['marker style'],
                   s=1.3*row['marker size']**2, zorder=100,
                   c=color, alpha=0.9, edgecolors='none')
        ax.scatter(row['inclination'], row['R0/pc']*FIX,
                   marker=row['marker style'],
                   s=(row['marker size'] - 1.2)**2, zorder=100,
                   c='w', alpha=1.0, edgecolors='none')

ax.legend(ncol=1, fontsize='small',
          title='Simulation\ncontact discontinuity\n(Meyer et al. 2017)',
          frameon=True, loc="upper left").get_title().set_size('small')
ax.set(
    xlim=[-5, 95],
    xticks=[0, 15, 30, 45, 60, 75, 90],
    ylim=[0., None],
    xlabel=r"Inclination, $|i|$",
    ylabel=r"Projected apex distance: $R_{0}'$, parsec",
)        
sns.despine(trim=True)
fig.tight_layout(pad=0.5)
fig.savefig(plotfile)
print(plotfile, end='')
