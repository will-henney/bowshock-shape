import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from astropy.table import Table
from dust_blorentz_ode import streamline

try:
    thB_degrees = float(sys.argv[1])
    LFAC = float(sys.argv[2])
    MACH_ALFVEN = float(sys.argv[3])
    try: 
        ALPHA_DRAG = float(sys.argv[4])
    except:
        ALPHA_DRAG = 0.0
except:
    sys.exit(f"Usage: {sys.argv[0]} FIELD_ANGLE LORENTZ_FORCE MACH_ALFVEN [ALPHA_DRAG]")


suffix = f"b{int(thB_degrees):02d}-L{int(100*LFAC):04d}-Ma{int(10*MACH_ALFVEN):04d}"
if ALPHA_DRAG > 0.0:
    suffix += f"-alpha{int(100*ALPHA_DRAG):04d}"
figfile = sys.argv[0].replace('.py', f'-{suffix}.jpg')

sns.set_style('white')
sns.set_color_codes()
fig, ax = plt.subplots(figsize=(5, 5))
ny0, nz0 = 101, 51
y0grid = 0.001 + np.linspace(-5.0, 5.0, ny0)
z0grid = 0.001 + np.linspace(-2.5, 2.5, nz0)

xmin, xmax = [-4.99, 4.99]
ymin, ymax = [-4.99, 4.99]
zmin, zmax = [-4.99, 4.99]
xx, yy, zz, ww = [], [], [], []
xs, ys, zs = [], [], []
nt = 3001
tstop = 60
X0 = 20
v_turb_0 = 1.0/MACH_ALFVEN
for iy0, y0 in enumerate(y0grid):
    for iz0, z0 in enumerate(z0grid):
        s = streamline(X0=X0, Y0=y0, Z0=z0, thB=np.radians(thB_degrees),
                       tstop=tstop, n=nt,
                       V_TURB_0=v_turb_0, LFAC=LFAC, ALPHA_DRAG=ALPHA_DRAG)
        # ax.plot(s['x'], s['y'], color='k', lw=0.5)
        # Accumulate (x, y) points in a long list
        xx.extend(s['x'])
        yy.extend(s['y'])
        zz.extend(s['z'])
        ww.extend(1.0 / (s['z']**2 + s['x']**2 + s['y']**2)**0.5)
        # ax.plot(s['x'], s['y'], '.',
        #         mec='none', mfc='r', ms=3, alpha=0.02)
        if iy0 % 1 == 0 and iz0 in [nz0//2,]:
            # Save streamlines for selected impact parameters
            xs.append(s['x'])
            ys.append(s['y'])
            zs.append(s['z'])
# Plot a density histogram of all the (x, y) points we accumulated
H, xe, ye = np.histogram2d(xx, yy, bins=(100, 100), weights=ww,
                           range=[[xmin, xmax], [ymin, ymax]])
# Do another 3d one with uniform weights for just the density
Hd, [ze, ye, ze] = np.histogramdd(np.stack((xx, yy, zz), axis=1),
                                  bins=(100, 100, 100), 
                                  range=[[xmin, xmax], [ymin, ymax], [zmin, zmax]])
np.savez(figfile.replace('.jpg', ''), rho=Hd,
         x=0.5*(xe[1:]+xe[:-1]), y=0.5*(ye[1:]+ye[:-1]), z=0.5*(ze[1:]+ze[:-1]))

rho_m = np.median(H)
ax.imshow(H.T, origin='lower', extent=[xmin, xmax, ymin, ymax],
          vmin=0.0, vmax=10.0*rho_m, cmap='gray_r')
# Plot the streamlines that we saved earlier
x1s = [4, 2, 0, -2, -4]
# Bespoke collection of colors for the scatter shot grains at each x1 position
colors = [
    # First batch at x=4 is separated from rest
    # Pale yellow
    (0.8, 0.8, 0.3, 1.0),
    # Swing towards orange for next two at x=2 and x=0
    (0.9, 0.7, 0.2, 1.0),
    (1.0, 0.4, 0.1, 1.0),
    # Now we need a contrast - go more purple
    (0.7, 0.1, 0.4, 1.0),
    # And finally, more blue
    (0.2, 0.2, 0.5, 1.0)
]
for x, y, z in zip(xs, ys, zs):
    ax.plot(x, y, '-', color='w', lw=0.8, alpha=0.5)
    ax.plot(x, y, '-', color='k', lw=0.5)
    for x1, color in zip(x1s, colors):
        itime = int((X0 - x1)*nt/tstop)
        ax.plot(x[itime-10:itime+20:10], y[itime-10:itime+20:10], '.', ms=4.0+z[itime], color=color, zorder=20-x1)
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


ax.set(
    ylabel='$y/R_{0}$',
    ylim=[ymin, ymax],
    xlabel='$x/R_{0}$',
    xlim=[xmin, xmax])

sns.despine()
fig.tight_layout()
fig.savefig(figfile, dpi=600)
print(figfile, end='')
