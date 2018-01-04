import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from astropy.table import Table
from dust_bfield_ode import streamline

figfile = sys.argv[0].replace('.py', '.jpg')

sns.set_style('white')
sns.set_color_codes()
fig, ax = plt.subplots(figsize=(5, 5))
ny0 = 5*400 + 1
y0grid = 0.001 + np.linspace(-5.0, 5.0, ny0)
nth = 200
thm_grid = np.linspace(0.0, np.pi, nth)
dth = np.pi/nth

rm = 2.0/(1.0 + np.cos(thm_grid))
xlocus = rm*np.cos(thm_grid)
ylocus = rm*np.sin(thm_grid)
xmin, xmax = [-4.99, 4.99]
ymin, ymax = [-4.99, 4.99]
xx, yy, ww = [], [], []
xs, ys = [], []
Z0 = 0.0
thB_degrees = 10.0
for iy0, y0 in enumerate(y0grid):
    s = streamline(X0=20, Y0=y0, Z0=Z0, thB=np.radians(thB_degrees), tstop=60, n=3001)
    # ax.plot(s['x'], s['y'], color='k', lw=0.5)
    # Accumulate (x, y) points in a long list
    xx.extend(s['x'])
    yy.extend(s['y'])
    ww.extend(1.0 / (Z0**2 + s['x']**2 + s['y']**2)**0.5)
    # ax.plot(s['x'], s['y'], '.',
    #         mec='none', mfc='r', ms=3, alpha=0.02)
    if iy0 % 30 == 15:
        # Save streamlines for selected impact parameters
        xs.append(s['x'])
        ys.append(s['y'])
# Plot a density histogram of all the (x, y) points we accumulated
H, xe, ye = np.histogram2d(xx, yy, bins=(100, 100), weights=ww,
                           range=[[xmin, xmax], [ymin, ymax]])
rho_m = np.median(H)
ax.imshow(H.T, origin='lower', extent=[xmin, xmax, ymin, ymax],
          vmin=0.0, vmax=20.0*rho_m, cmap='gray_r')
# Plot the streamlines that we saved earlier
for x, y in zip(xs, ys):
    ax.plot(x, y, '-', color='w', lw=0.8, alpha=0.5)
    ax.plot(x, y, '-', color='k', lw=0.5)
ax.plot(xlocus, ylocus, ':', color='y', alpha=0.7, lw=2)
ax.plot(xlocus, -ylocus, ':', color='y', alpha=0.7, lw=2)
cthB = np.cos(np.radians(thB_degrees))
sthB = np.sin(np.radians(thB_degrees))
for xx in np.linspace(1.5*xmin, 1.5*xmax, 15):
    yy1, yy2 = 1.5*ymin, 1.5*ymax
    x1 = -xx*sthB + yy1*cthB
    x2 = -xx*sthB + yy2*cthB
    y1 = xx*cthB + yy1*sthB
    y2 = xx*cthB + yy2*sthB
    ax.plot([x1, x2], [y1, y2], lw=2, alpha=0.5, color='c')
ax.axvline(0.0, ls='--', color='0.5', lw=0.5)
ax.axhline(0.0, ls='--', color='0.5', lw=0.5)
ax.plot([0.0], [0.0], '+', color='k')
ax.plot([0.0], [0.0], '.', color='k')
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

tabfilename = sys.argv[0].replace('.py', f'.tab')
Table({'theta': thm_grid, 'R': rrm_grid}).write(tabfilename, format='ascii.tab')

ax.set(
    ylabel='$y/R_{0}$',
    ylim=[ymin, ymax],
    xlabel='$x/R_{0}$',
    xlim=[xmin, xmax])

sns.despine()
fig.tight_layout()
fig.savefig(figfile, dpi=600)
print(figfile, end='')
