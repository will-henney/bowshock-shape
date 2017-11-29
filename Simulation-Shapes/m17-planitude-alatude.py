T=[["Source", "R0/100", "Rc/R0", "R90/R0", "Rm90/R0"], ["M17-MHD2040-AllB7", "0.48", "3.57", "1.95", "1.92"], ["M17-MHD2040-AllB7-Halpha-i00", "0.62", "4.78", "2.21", "2.21"], ["M17-MHD2040-AllB7-60mic-i30", "0.20", "2.14", "1.74", "1.73"], ["M17-MHD2040-AllB7-60mic-i45", "0.23", "1.70", "1.59", "1.61"], ["M17-MHD2040-AllB7-60mic-i60", "0.27", "1.51", "1.47", "1.46"], ["M17-HD2040", "0.84", "1.78", "1.75", "1.75"], ["M17-HD2040-Halpha-i00-CD", "0.99", "1.96", "2.09", "2.08"], ["M17-HD2040-Halpha-i00-BS", "1.31", "2.15", "1.95", "1.96"], ["M17-HD2040-60mic-i30", "0.48", "1.89", "1.80", "1.79"], ["M17-HD2040-60mic-i45", "0.56", "1.87", "1.90", "1.59"], ["M17-HD2040-60mic-i60", "0.70", "1.66", "1.55", "1.59"]]
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
    else:
        return 5

table['marker style'] = [select_marker_style(s) for s in table['Source']]
table['marker size'] = [select_marker_size(s) for s in table['Source']]

sns.set_style('ticks')
fig, ax = plt.subplots(figsize=(4, 4))

Rc_grid = np.linspace(0.9, 10.0, 2000)
R90_T0_grid = np.sqrt(2*Rc_grid)
R90_T1_grid = np.sqrt(2*Rc_grid - 1.0)
R90_T1_grid[~np.isfinite(R90_T1_grid)] = 0.0 

ax.fill_between(Rc_grid, R90_T1_grid, R90_T0_grid, color='k', alpha=0.2)
ax.fill_between(Rc_grid, R90_T0_grid, 0.5, color='k', alpha=0.1)
ax.plot(Rc_grid, R90_T0_grid, c='k', lw=0.5)
ax.axhline(1.0, lw=0.5, alpha=0.5, color='k', zorder=-1)
ax.axvline(1.0, lw=0.5, alpha=0.5, color='k', zorder=-1)
ax.plot([0.9, 10.0], [0.9, 10.0], lw=0.5, alpha=0.5, color='k', zorder=-1)

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



bp.N_NEIGHBORHOOD = 200
bp.DEGREE_POLY_NEIGHBORHOOD = 1
bp.SCALE_NEIGHBORHOOD = (60./180.) # => +/-60 deg at i=0
bp.DEGREE_POLY_NEIGHBORHOOD_90 = 2
bp.SCALE_NEIGHBORHOOD_90 = 0.01

models = ["M17-HD2040", "M17-MHD2040-AllB7"][::-1]
labels = ["HD", "MHD"][::-1]

colors = sns.color_palette(n_colors=len(models))[::-1]
for model, label, color in zip(models, labels, colors):
    # First do the "theoretical" tracks
    shape = Simulation(name=model, force_open=True, cheby_degree=12)
    th_inf = bp.theta_infinity(shape)
    inc = np.linspace(0.0, th_inf - np.pi/2, 50)
    tab = bow_diagnostic.parameter_table(inc, shape)
    Rc, R90 = tab['tilde R_c prime'], tab['tilde R_90 prime']
    ax.plot(Rc, R90, '-', c=color, label=label, lw=2.5, alpha=0.5)
    # Put a dot at the i=0 case
    ax.plot(Rc[0:1], R90[0:1], 'o', mec='none', c=color, label="_nolabel_", alpha=0.5)

    # Repeat for the closed shape, but with a thin line
    shape = Simulation(name=model, force_open=False, cheby_degree=12)
    th_inf = bp.theta_infinity(shape)
    inc = np.linspace(0.0, th_inf - np.pi/2, 50)
    tab = bow_diagnostic.parameter_table(inc, shape)
    Rcc, R90c = tab['tilde R_c prime'], tab['tilde R_90 prime']
    ax.plot(Rcc, R90c, '-', c=color, label="_nolabel_", lw=0.5, alpha=1.0)


    sini = (0.5 + np.arange(20))/20
    inc_e = np.arcsin(sini)
    inc_e = inc_e[inc_e < th_inf - np.pi/2]
    # Interpolate to get the even probability points
    Rc_e = interp1d(inc, Rc)(inc_e)
    R90_e = interp1d(inc, R90)(inc_e)

    ax.scatter(Rc_e, R90_e, marker='|', s=3**2,
               linewidths=0.1, edgecolors='none',
               c=color, alpha=0.5, label="_nolabel_")

    # Label i=0, 30, 45, 60 along the MHD track
    interp_Rc = interp1d(inc, Rc, bounds_error=False, fill_value="extrapolate")
    interp_R90 = interp1d(inc, R90, bounds_error=False, fill_value="extrapolate")
    if "MHD" in model:
        annotate_data = [
            [0, (-20, -25)], [30, (15, -20)],
            [45, (5, -25)], [60, (-10, -30)]]
    else:
        annotate_data = [
            [0, (-20, 25)], [30, (15, 15)],
            [60, (-20, 25)]]

    for inclination, xytext in annotate_data:
        Rc_i = interp_Rc(np.radians(inclination))
        R90_i = interp_R90(np.radians(inclination))
        ax.annotate(fr"$i = {inclination}^\circ$",
                    (Rc_i, R90_i),
                    xytext=xytext, textcoords='offset points',
                    arrowprops=dict(arrowstyle="->", color=color,
                                    connectionstyle="arc3,rad=.2"),       
                    color=color, fontsize="x-small")


    mask = [s.startswith(model) for s in table['Source']]
    data = table[mask]

    for row in data:
        ax.scatter(row['Rc/R0'], row['R90'],
                   marker=row['marker style'],
                   s=1.3*row['marker size']**2, zorder=100,
                   c=color, alpha=0.9, edgecolors='none')
        ax.scatter(row['Rc/R0'], row['R90'],
                   marker=row['marker style'],
                   s=0.5*row['marker size']**2, zorder=100,
                   c='w', alpha=1.0, edgecolors='none')


ax.legend(ncol=1, fontsize='small',
          title='Simulations\n(Meyer et al. 2017)',
          frameon=True, loc="upper left").get_title().set_size('small')
ax.set(
    xlim=[0.93, 6.1],
    ylim=[0.93, 6.1],
    yscale='log',
    xscale='log',
    #yticks=range(6),
    xlabel=r"Projected planitude: $\Pi'$",
    ylabel=r"Projected alatude: $\Lambda'$",
)        
ax.xaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter('%d'))
ax.yaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter('%d'))
ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%d'))
ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%d'))
sns.despine()
fig.tight_layout(pad=0.5)
fig.savefig(plotfile)
print(plotfile, end='')
