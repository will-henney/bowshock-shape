import sys
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from astropy.table import Table
from dust_blorentz_ode import streamline

try:
    thB_degrees = float(sys.argv[1])
    LFAC = float(sys.argv[2])
except:
    sys.exit(f"Usage: {sys.argv[0]} FIELD_ANGLE LORENTZ_FORCE")


suffix = f"b{int(thB_degrees):02d}-L{int(100*LFAC):04d}"
figfile = sys.argv[0].replace('.py', f'-{suffix}.jpg')

sns.set_style('white')
sns.set_color_codes()
fig, ax = plt.subplots(subplot_kw=dict(projection='3d'), figsize=(6, 3.2))


ny0, nz0 = 31, 1
y0grid = 0.000 + np.linspace(-2.5, 2.5, ny0)
z0grid = [0.0]

xmin, xmax = [-6, 2]
zmin, zmax = [-2, 2]
ymin, ymax = [-3.5, 3.5]
#nt = 3001
nt = 3001
tstop = 60
X0 = 20

zcolors = ['k']
ycolors = sns.color_palette('magma_r', n_colors=ny0)

for iz0, z0 in enumerate(z0grid):
    for iy0, y0 in enumerate(y0grid[::-1]):
        s = streamline(X0=X0, Y0=y0, Z0=z0, thB=np.radians(thB_degrees),
                       tstop=tstop, n=nt, LFAC=LFAC)
        x, y, z = s['x'], s['y'], s['z']
        # implement clipping by hand
        m = (x >= xmin) & (x <= xmax)
        m = m & (y >= ymin) & (y <= ymax)
        m = m & (z >= zmin) & (z <= zmax)
        x[~m] = np.nan
        if iy0 == ny0//2:
            ax.plot([0, 0], [0, 0], [zmin, zmax], '--', color='k', lw=0.3)
            ax.plot([xmin, xmax], [0, 0], [0, 0], '--', color='k', lw=0.3)
            ax.plot([0, 0], [ymin, ymax], [0, 0], '--', color='k', lw=0.3)
            star = ax.plot([0], [0], [0], 'o', color='r')
        ax.plot(x, y, z, '-', color='w', lw=1.0)
        ax.plot(x, y, z, '-', color=ycolors[iy0], lw=0.7)

cthB = np.cos(np.radians(thB_degrees))
sthB = np.sin(np.radians(thB_degrees))
for xx in np.linspace(-15, 15, 31):
    yy1, yy2 = -15, 15
    x1 = -xx*sthB + yy1*cthB
    x2 = -xx*sthB + yy2*cthB
    y1 = xx*cthB + yy1*sthB
    y2 = xx*cthB + yy2*sthB
    x = np.linspace(x1, x2, 200)
    y = np.linspace(y1, y2, 200)
    m = (x >= xmin) & (y >= ymin) & (x <= xmax) & (y <= ymax)
    x[~m] = np.nan
    ax.plot(x, y, zs=zmin, zdir='z', lw=2, alpha=0.5, color='c')

ax.auto_scale_xyz([xmin, xmax], [ymin, ymax], [zmin, zmax])
ax.set(
    xlabel='$x/R_{0}$',
    ylabel='$y/R_{0}$',
    zlabel='$z/R_{0}$',
    zticks=[-2, -1, 0, 1, 2],
    xlim=[xmin, xmax],
    ylim=[ymin, ymax],
    zlim=[zmin, zmax],
)
ax.azim = -55
ax.elev = 20
fig.text(0.1, 0.8, fr"$\theta_B = {thB_degrees:.0f}^\circ \quad\quad F_B / F_\mathrm{{rad}} = {LFAC:.1f}$")
fig.tight_layout(rect=[0, 0, 0.95, 1.0])
fig.savefig(figfile, dpi=600)
print(figfile, end='')
