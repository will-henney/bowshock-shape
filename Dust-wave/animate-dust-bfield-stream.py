import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import seaborn as sns
from astropy.table import Table
from dust_bfield_ode import streamline

try:
    thB_degrees = float(sys.argv[1])
    Z0 = float(sys.argv[2])
except:
    sys.exit(f"Usage: {sys.argv[0]} FIELD_ANGLE Z0")


suffix = f"b{int(thB_degrees):02d}-z{int(100*Z0):04d}"
fileroot = sys.argv[0].replace('.py', f'-{suffix}')

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
tstop = 60
nt = 3001
X0 = 20
for iy0, y0 in enumerate(y0grid):
    s = streamline(X0=X0, Y0=y0, Z0=Z0, thB=np.radians(thB_degrees), tstop=tstop, n=nt)
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

# Artist to use for animation of streamlines
point_artists = [ax.plot([], [], '.', color='r')[0] for _ in xs]
nB = 30
line_artists = [ax.plot([], [], lw=2, alpha=0.5, color='c')[0] for _ in range(nB)]
xxBs = np.linspace(2*xmin, 2*xmax, nB)

def animate_points_and_lines(itime):
    """Update time along all plotted streamline points and lines"""
    # Update the lines
    for xx, artist in zip(xxBs, line_artists):
        yy1, yy2 = 3*ymin, 3*ymax
        x1 = -xx*sthB + yy1*cthB  + X0 - tstop*itime/nt
        x2 = -xx*sthB + yy2*cthB  + X0 - tstop*itime/nt
        y1 = xx*cthB + yy1*sthB
        y2 = xx*cthB + yy2*sthB
        artist.set_data([x1, x2], [y1, y2])
    # Update the points
    for x, y, artist in zip(xs, ys, point_artists):
        artist.set_data(x[itime-10:itime+20:10], y[itime-10:itime+20:10])
    return point_artists + line_artists

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
ax.axvline(0.0, ls='--', color='0.5', lw=0.5)
ax.axhline(0.0, ls='--', color='0.5', lw=0.5)
ax.plot([0.0], [0.0], '+', color='k')
ax.plot([0.0], [0.0], '.', color='k')
ax.set_aspect('equal', adjustable='box-forced')

ax.set(
    ylabel='$y/R_{0}$',
    ylim=[ymin, ymax],
    xlabel='$x/R_{0}$',
    xlim=[xmin, xmax])

sns.despine()
fig.tight_layout()

anim = animation.FuncAnimation(fig, animate_points_and_lines,
                               frames=range(700, 1501, 2),
                               blit=True)
moviefile = fileroot + '.mp4'
anim.save(moviefile, writer='ffmpeg', fps=60, dpi=200)
print(moviefile, end='')
