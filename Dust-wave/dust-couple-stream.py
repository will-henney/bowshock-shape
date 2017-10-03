import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from dust_couple_ode import streamline

figfile = sys.argv[0].replace('.py', '.jpg')

sns.set_style('white')
sns.set_color_codes()
fig, axes = plt.subplots(2, 2, figsize=(6, 4))
alphas = [1.0/4.0, 1.0/2.0, 1.0, 2.0]
nb = 501
bgrid = 0.005 + np.linspace(0.0, 5.0, 501)
ibspecial = [5, 20, 50, 100, 200, 300]
thm_grid = np.linspace(0.0, np.pi, 200)
rm = 2.0/(1.0 + np.cos(thm_grid))
xlocus = rm*np.cos(thm_grid)
ylocus = rm*np.sin(thm_grid)
xmin, xmax = [-3.999, 3.999]
ymin, ymax = [0.0, 4.999]
for alpha, ax in zip(alphas, axes.flat):
    xx, yy = [], []
    xs, ys = [], []
    for ib, b in enumerate(bgrid):
        s = streamline(y0=b, alpha=alpha, n=5001)
        # ax.plot(s['x'], s['y'], color='k', lw=0.5)
        # Accumulate (x, y) points in a long list
        xx.extend(s['x'])
        yy.extend(s['y'])
        # ax.plot(s['x'], s['y'], '.',
        #         mec='none', mfc='r', ms=3, alpha=0.02)
        if ib in ibspecial:
            # Save streamlines for selected impact parameters
            xs.append(s['x'])
            ys.append(s['y'])
    # Plot a density histogram of all the (x, y) points we accumulated
    ax.hist2d(xx, yy, bins=(80*2, 50*2),
              range=[[xmin, xmax], [ymin, ymax]])
    # Plot the streamlines that we saved earlier
    for x, y in zip(xs, ys):
        ax.plot(x, y, '-', color='y', lw=0.5)
    ax.plot(xlocus, ylocus, ':', color='w', alpha=0.5, lw=2)
    ax.axvline(0.0, ls='--', color='w', lw=0.5)
    ax.text(1.0, 4.0, 
            fr"$\alpha_\mathrm{{drag}} = {alpha:.2f}$",
            color='y')
    ax.set(xlabel='$x$', ylabel='$y$',
           xlim=[xmin, xmax], ylim=[ymin, ymax])
    ax.set_aspect('equal', adjustable='box')


sns.despine()
fig.tight_layout()
fig.savefig(figfile, dpi=300)
print(figfile, end='')
