import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
from regress_utils import plot_regression

sns.set_color_codes('dark')

figfile = 'mipsgal-h0-thb.pdf'

combo_file = 'mipsgal-arcfit.tab'
tab = Table.read(combo_file, format='ascii.tab')

fig, ax = plt.subplots(figsize=(6, 6))

Q = tab['R0_fit']/tab['R0_sigma']

# Make quality classes
m5 = tab['Rating'] ==  5
m4 = tab['Rating'] ==  4
m3 = tab['Rating'] ==  3
m2 = tab['Rating'] ==  2
m1 = tab['Rating'] ==  1

masks = (m1 | m2), m3, m4, m5
alphas = 0.05, 0.3, 0.5, 0.8
labels = '1- or 2-star', '3-star', '4-star', '5-star'
colors = 'k', 'c', 'r', 'b'
sizes = 3, 5, 7, 10

range_ = [0.25, 8]
ticks = [0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]

# Add some columns
tab["th_b"] = 0.5*(tab["th_p"] - tab["th_m"])
tab["h"] = tab["H_g"]/tab["R0_g"]

mfinite = np.isfinite(tab["h"]) & np.isfinite(tab["th_b"])

for m, alpha, c, ms, label in zip(masks, alphas, colors, sizes, labels):
    ax.plot('h', "th_b", 'o', data=tab[m],
                lw=1, c=c, ms=ms, alpha=alpha,
                label=f'{label} ($N = {m.sum()}$)')


for source in tab[m4]:
    ax.text(source['h'], source['th_b'], f'{source["Seq"]}',
            fontsize=4, color='white', ha='center', va='center')
for source in tab[m5]:
    ax.text(source['h'], source['th_b'], f'{source["Seq"]}',
            fontsize=5, color='orange', ha='center', va='center')

#mgood = m3 | m4 | m5
mgood = (m5 | m4 ) & mfinite
# mgood = (m5 | m4 | m3 | m2 | m1) & mfinite
plot_regression(ax, tab[mgood]["h"], tab[mgood]["th_b"],
                logx=True, xlim=range_, debug=True, pos="top right")

ax.legend()
ax.set(
    xlim=range_, ylim=[0, 135],
    xlabel='Dimensionless shell thickness, $h$',
    ylabel='Brightness half width half max, degrees',
    xscale='log', yscale='linear')
ax.xaxis.set_major_locator(mticker.FixedLocator(ticks))
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%0g'))
ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%0g'))

sns.despine()
fig.savefig(figfile)
