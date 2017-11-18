import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from astropy.table import Table
from dust_couple_ode import streamline

figfile = sys.argv[0].replace('.py', '.jpg')

sns.set_style('ticks')
sns.set_color_codes()
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(6, 4))
alphas = [1.0, 1.0, 4.0, 4.0]
mus = [0.05, 0.2, 0.2, 0.8]
nb = 25*200 + 1
bgrid = 0.001 + np.linspace(0.0, 5.0, nb)
ibspecial = 25*np.array([1, 4, 10, 20, 40, 60])
nth = 200
thm_grid = np.linspace(0.0, np.pi, nth)
dth = np.pi/nth


xmin, xmax = [-4.1, 4.1]
ymin, ymax = [0.0, 5.1]
for alpha, mu, ax in zip(alphas, mus, axes.flat):
    xx, yy, ww = [], [], []
    xs, ys = [], []

    # zoom in on the alpha > 1 models since they get small
    zoom = alpha if alpha > 1.0 else 1.0

    # Launch grains on a uniform grid of th1
    # Make sure it fills the plot
    th1max = np.arctan2(ymax/zoom, 1.0/mu - xmax/zoom)
    th1grid = 0.001*mu + np.linspace(0.0, th1max, nb)
    bgrid = np.sin(th1grid)/mu

    # Hyperbola solution for drag-free case, but scaling mu by alpha
    ecc = 1.0 / (1.0 - 2*mu/alpha)
    # And scale radius by alpha too
    rm = (1.0 + ecc)/(1.0 + ecc*np.cos(thm_grid))/alpha
    rm[rm < 0.0] = np.nan
    xlocus = rm*np.cos(thm_grid)
    ylocus = rm*np.sin(thm_grid)


    for ib, (th1, b) in enumerate(zip(th1grid, bgrid)):
        # Start from a circle just outside the plot window
        Rlaunch = 1/mu - xmax/zoom
        assert Rlaunch > 0.0
        X0 = 1./mu - Rlaunch*np.cos(th1)
        Y0 = Rlaunch*np.sin(th1)
        s = streamline(X0=X0, Y0=Y0, tstop=30, alpha=alpha, mu=mu, n=30001)
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
                               range=[[xmin/zoom, xmax/zoom], [ymin/zoom, ymax/zoom]])
    rho_m = np.median(H[H != 0.0])
    rho_m = H[-1, -1]
    ax.imshow(H.T, origin='lower', extent=[xmin, xmax, ymin, ymax],
              vmin=0.0, vmax=2.0*rho_m, cmap='gray_r')
    # Plot the streamlines that we saved earlier
    for x, y in zip(xs, ys):
        ax.plot(x*zoom, y*zoom, '-', color='w', lw=0.8, alpha=0.5)
        ax.plot(x*zoom, y*zoom, '-', color='k', lw=0.5)
    ax.plot(xlocus*zoom, ylocus*zoom, ':', color='w', alpha=0.5, lw=2)
    ax.axvline(0.0, ls='--', color='w', lw=0.5)
    label = fr"$\alpha_\mathrm{{drag}} = {alpha:.1f}$"
    label += '\n' + fr"$\mu = {mu:.2f}$"
    ax.text(1.0, 4.0, label, color='k')
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

    suffix = f'-alpha{int(100*alpha):03d}'
    suffix += f'-mu{int(100*mu):03d}'
    tabfilename = sys.argv[0].replace('.py', suffix + '.tab')
    Table({'theta': thm_grid, 'R': rrm_grid}).write(tabfilename, format='ascii.tab', overwrite=True)

for ax in axes[:, 0]:
    ax.set(
        ylabel=r'$\alpha_\mathrm{{drag}} \,Y$',
        ylim=[ymin, ymax],
        yticks=range(5),
    )
for ax in axes[-1, :]:
    ax.set(
        xlabel=r'$\alpha_\mathrm{{drag}} \,X$',
        xlim=[xmin, xmax],
        xticks=range(-4,5),
    )

sns.despine()
fig.tight_layout()
fig.savefig(figfile, dpi=600)
print(figfile, end='')
