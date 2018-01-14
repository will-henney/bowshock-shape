import sys
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Rectangle
import seaborn as sns
from astropy.table import Table
from dust_blorentz_ode import streamline

try:
    thB_degrees = float(sys.argv[1])
    LFAC = float(sys.argv[2])
    MACH_ALFVEN = float(sys.argv[3])
    YPLANE = float(sys.argv[4])
except:
    sys.exit(f"Usage: {sys.argv[0]} FIELD_ANGLE LORENTZ_FORCE MACH_ALFVEN YPLANE")


suffix = f"b{int(thB_degrees):02d}-L{int(100*LFAC):04d}-Ma{int(10*MACH_ALFVEN):04d}-Y{int(100*YPLANE):+04d}"
figfile = sys.argv[0].replace('.py', f'-{suffix}.jpg')

sns.set_style('white')
sns.set_color_codes()
fig, ax = plt.subplots(subplot_kw=dict(projection='3d'), figsize=(6, 3.3))


ny0, nz0 = 31, 1
y0grid = 0.000 + np.linspace(-2.0, 2.0, 1.5*ny0)
z0grid = [0.0]

xmin, xmax = [-7, 3]
zmin, zmax = [-2, 2]
ymin, ymax = [-4, 4]

# show the z=0 plane
p = Rectangle((xmin, zmin), xmax - xmin, zmax - zmin,
              edgecolor='none', facecolor='g', alpha=0.15)
ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=YPLANE, zdir="y")

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

#nt = 3001
nt = 3001
tstop = 60
X0 = 10
# Time index corresponding to stellar passage, approx
it0 = int(nt*X0/tstop)

zcolors = ['k']

# Random velocity of grains - assume at Alfven speed, at high pitch angle to field
v_turb_0 = 1.0/MACH_ALFVEN
mu_p_max = 0.1

# Find all the trajectories
trajectories = []
for y0 in y0grid:
    # Pitch angle to field
    c_p = np.random.uniform(0.0, mu_p_max)
    s_p = np.sqrt(1.0 - c_p**2) 
    # Azimuth around field
    phi_B = np.random.uniform(0.0, 2*np.pi)
    spB, cpB = np.sin(phi_B), np.cos(phi_B) 
    # Magnitude and sign
    v_turb = np.random.normal(loc=0.0, scale=v_turb_0)
    # Components in frame of B-field
    vBx = v_turb*c_p
    vBy = v_turb*s_p*cpB
    vBz = v_turb*s_p*spB
    # Now transform to flow frame
    cthB, sthB = np.cos(np.radians(thB_degrees)), np.sin(np.radians(thB_degrees))
    u0 = vBx*cthB - vBy*sthB - 1.0
    v0 = vBx*sthB + vBy*cthB
    w0 = vBz

    s = streamline(X0=X0, Y0=YPLANE, Z0=y0, U0=u0, V0=v0, W0=w0,
                   thB=np.radians(thB_degrees),
                   tstop=tstop, n=nt, LFAC=LFAC)
    if ymin <= s['y'][it0] <= ymax:
        trajectories.append(s)

# Sort according to z coordinate near the star
trajectories.sort(key=lambda s: s['z'][it0], reverse=False)

ycolors = sns.color_palette('magma_r', n_colors=len(trajectories))

star_done = False
for iy0, s in enumerate(trajectories):
    x, y, z = s['x'], s['y'], s['z']
    # implement clipping by hand
    m = (x >= xmin) & (x <= xmax)
    m = m & (y >= ymin) & (y <= ymax)
    m = m & (z >= zmin) & (z <= zmax)
    x[~m] = np.nan
    # Draw the star before we get to the first positive z value
    if s['z'][it0] > 0.0 and not star_done:
        ax.plot([0, 0], [0, 0], [zmin, zmax], '--', color='k', lw=0.3)
        ax.plot([xmin, xmax], [0, 0], [0, 0], '--', color='k', lw=0.3)
        ax.plot([0, 0], [ymin, ymax], [0, 0], '--', color='k', lw=0.3)
        star = ax.plot([0], [0], [0], 'o', color='r')
        star_done = True
    ax.plot(x, y, z, '-', color='w', lw=1.0)
    ax.plot(x, y, z, '-', color=ycolors[iy0], lw=0.7)


ax.auto_scale_xyz([xmin, xmax], [ymin, ymax], [zmin, zmax])
ax.set(
    xlabel='$x/R_{0}$',
    ylabel='$y/R_{0}$',
    zlabel='$z/R_{0}$',
    xticks=range(xmin, xmax+1),
    yticks=range(ymin, ymax+1),
    zticks=[-2, -1, 0, 1, 2],
    xlim=[xmin, xmax],
    ylim=[ymin, ymax],
    zlim=[zmin, zmax],
    )
ax.azim = +15
ax.elev = 30
text = "$" + r" \quad\quad ".join([
    fr"\theta_B = {thB_degrees:.0f}^\circ",
    fr"F_B / F_\mathrm{{rad}} = {LFAC:.1f}",
    fr"\mathcal{{M}}_A = {MACH_ALFVEN:.1f}",
    fr"y_0 = {YPLANE:.4f}"]) + "$"
fig.text(0.1, 0.9, text)
fig.tight_layout(rect=[0, 0, 1.05, 1.0])
fig.savefig(figfile, dpi=600)
print(figfile, end='')
