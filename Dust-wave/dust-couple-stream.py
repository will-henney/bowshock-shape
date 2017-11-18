import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from astropy.table import Table
from dust_couple_ode import streamline

figfile = sys.argv[0].replace('.py', '.jpg')

sns.set_style('white')
sns.set_color_codes()
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(6, 4))
alphas = [1.0/4.0, 1.0/2.0, 1.0, 2.0]
nb = 25*200 + 1
bgrid = 0.001 + np.linspace(0.0, 5.0, nb)
ibspecial = 25*np.array([1, 4, 10, 20, 40, 60])
nth = 200
thm_grid = np.linspace(0.0, np.pi, nth)
dth = np.pi/nth

rm = 2.0/(1.0 + np.cos(thm_grid))
xlocus = rm*np.cos(thm_grid)
ylocus = rm*np.sin(thm_grid)
xmin, xmax = [-3.99, 3.99]
ymin, ymax = [0.0, 4.99]
for alpha, ax in zip(alphas, axes.flat):
    xx, yy, ww = [], [], []
    xs, ys = [], []
    for ib, b in enumerate(bgrid):
        s = streamline(X0=5, Y0=b, tstop=30, alpha=alpha, n=30001)
        # ax.plot(s['x'], s['y'], color='k', lw=0.5)
        # Accumulate (x, y) points in a long list
        xx.extend(s['x'])
        yy.extend(s['y'])
        # Weights proportional to b/r
        ww.extend(s['b']/s['y'])
        # ax.plot(s['x'], s['y'], '.',
        #         mec='none', mfc='r', ms=3, alpha=0.02)
        if ib in ibspecial:
            # Save streamlines for selected impact parameters
            xs.append(s['x'])
            ys.append(s['y'])
    # Plot a density histogram of all the (x, y) points we accumulated
    H, xe, ye = np.histogram2d(xx, yy, bins=(80/1, 50/1), weights=ww,
                               range=[[xmin, xmax], [ymin, ymax]])
    rho_m = np.median(H)
    ax.imshow(H.T, origin='lower', extent=[xmin, xmax, ymin, ymax],
              vmin=0.0, vmax=2.0*rho_m, cmap='gray_r')
    # Plot the streamlines that we saved earlier
    for x, y in zip(xs, ys):
        ax.plot(x, y, '-', color='w', lw=0.8, alpha=0.5)
        ax.plot(x, y, '-', color='k', lw=0.5)
    ax.plot(xlocus, ylocus, ':', color='w', alpha=0.5, lw=2)
    ax.axvline(0.0, ls='--', color='w', lw=0.5)
    ax.text(1.0, 4.0, 
            fr"$\alpha_\mathrm{{drag}} = {alpha:.2f}$",
            color='k')
    ax.set_aspect('equal', adjustable='box-forced')

    # Save the minimum radius as a function of theta
    rr = np.hypot(xx, yy)
    theta = np.arctan2(yy, xx)
    rrm_grid = np.empty_like(thm_grid)
    for j, th0 in enumerate(thm_grid):
        # Mask to select points with theta between th0 -> th0 + dth
        m = np.abs(theta - (th0 + 0.5*dth)) <= 0.5*dth
        try:
            rrm_grid[j] = rr[m].min()
        except:
            # Sometimes mask may be empty
            rrm_grid[j] = np.nan

    tabfilename = sys.argv[0].replace('.py', f'-alpha{int(100*alpha):03d}.tab')
    Table({'theta': thm_grid, 'R': rrm_grid}).write(tabfilename, format='ascii.tab')

for ax in axes[:, 0]:
    ax.set(ylabel='$y/R_{0}$', ylim=[ymin, ymax])
for ax in axes[-1, :]:
    ax.set(xlabel='$x/R_{0}$', xlim=[xmin, xmax])

sns.despine()
fig.tight_layout()
fig.savefig(figfile, dpi=600)
print(figfile, end='')
