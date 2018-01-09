import sys
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import seaborn as sns

try:
    ANGLE = int(sys.argv[1])
except:
    sys.exit(f"Usage: {sys.argv[0]} ANGLE")


figfile = sys.argv[0].replace('.py', f'-{ANGLE:02d}.pdf')
fig, axes = plt.subplots(3, 4, sharex=True, sharey=True, figsize=(8, 6.1))

# Unit vector along B-field in the body-frame basis
cth = np.cos(np.radians(ANGLE))
sth = np.sin(np.radians(ANGLE))
bvec = np.array([cth, sth, 0.0])

# Unit vector along flow direction in the body-frame basis
uvec = np.array([-1.0, 0.0, 0.0])

def circle(x, y, radius=0.15):
    """Copied from 'Anatomy of a figure' example in matplotlib docs"""
    from matplotlib.patches import Circle
    from matplotlib.patheffects import withStroke
    circle = Circle((x, y), radius, clip_on=False, linewidth=0.3, ls=':',
                    edgecolor='w', facecolor=(1, 1, 1, .125),
                    path_effects=[withStroke(linewidth=0.8, foreground='0.5')])
    ax.add_artist(circle)


xmin, xmax, ymin, ymax = -5.0, 5.0, -5.0, 5.0
for j, inc in enumerate([7, 22, 39, 61]):
    for i, azi in enumerate([15, 45, 75]):
        fitsfile = f"dust-bfield-b{ANGLE:02d}-cube-F-{inc:03d}-{azi:03d}.fits"
        hdu, = fits.open(fitsfile)

        ax = axes[i, j]
        ax.imshow(np.sqrt(hdu.data),
                  origin='lower', extent=[xmin, xmax, ymin, ymax])
        thBtext = fr"$\theta_B = {ANGLE}^\circ$"
        inctext = fr"$|i| = {inc}^\circ$"
        azitext = fr"$|\phi| = {azi}^\circ$"
        ax.text(-4.2, -4.2, inctext, color="w")
        ax.text(4.2, -4.2, azitext, color="w", ha='right')
        ax.text(4.2, 4.2, thBtext, color="w", ha='right', va='top')
        ax.plot([0], [0], '+', color='orange')

        ci = np.cos(np.radians(inc))
        si = np.sin(np.radians(inc))
        ca = np.cos(np.radians(azi))
        sa = np.sin(np.radians(azi))

        # Observer-frame unit vectors in the body-frame basis
        iprime = np.array([ci, -si*sa, -si*ca])
        jprime = np.array([0.0, ca, -sa])
        kprime = np.array([si, ci*sa, ci*ca])

        # Observer frame components of B-field
        bx = np.dot(bvec, iprime)
        by = np.dot(bvec, jprime)
        bz = np.dot(bvec, kprime)

        bcolor = (0.2, 0.9, 1.0)
        bsky = np.hypot(bx, by)
        btheta = np.arctan2(by, bx)
        x0, y0, scale = -3.0, 3.0, 0.5
        circle(x0, y0, scale)
        ax.plot(
            [x0 - 1.1*bx*scale/bsky, x0 + 1.1*bx*scale/bsky],
            [y0 - 1.1*by*scale/bsky, y0 + 1.1*by*scale/bsky],
            lw=0.3, color='w', ls=':',
        )
        ax.arrow(x0 - bx*scale, y0 - by*scale,
                 2*bx*scale, 2*by*scale,
                 width=0.001, head_width=0.1,
                 head_length=0.1,
                 length_includes_head=True,
                 color=bcolor,
                 zorder=10,
        )
        btext = ax.annotate(r'$\vec{B}$', (x0, y0),
                            ha='center', color=bcolor,
                            xytext=(0,10), textcoords='offset points')
        btext.set_path_effects([
            path_effects.Stroke(linewidth=1, foreground=(0.0, 0.0, 0.0, 0.1)),
            path_effects.Normal()
        ])

        # Observer frame components of Velocity
        ux = np.dot(uvec, iprime)
        uy = np.dot(uvec, jprime)
        uz = np.dot(uvec, kprime)

        ucolor = (1.0, 0.8, 0.2)
        usky = np.hypot(ux, uy)
        x0, y0, scale = 3.0, 0.0, 0.5
        circle(x0, y0, scale)
        ax.plot(
            [x0 - 1.1*ux*scale/usky, x0 + 1.1*ux*scale/usky],
            [y0 - 1.1*uy*scale/usky, y0 + 1.1*uy*scale/usky],
            lw=0.3, color='w', ls=':',
        )
        ax.arrow(x0 - ux*scale, y0 - uy*scale,
                 2*ux*scale, 2*uy*scale,
                 width=0.001, head_width=0.1,
                 head_length=0.1,
                 length_includes_head=True,
                 color=ucolor,
                 zorder=10,
        )
        utext = ax.annotate(r'$\vec{v}_{\infty}$', (x0, y0),
                            ha='center', color=ucolor,
                            xytext=(0,10), textcoords='offset points')
        utext.set_path_effects([
        path_effects.Stroke(linewidth=1, foreground=(0.0, 0.0, 0.0, 0.1)),
            path_effects.Normal()
        ])



axes[-1, 0].set(
    xlabel="$x' / R_{_{**}}$",
    ylabel="$y'/ R_{_{**}}$",
    xlim=[-4.5, 4.5],
    ylim=[-4.5, 4.5],
    xticks=[-4, -2, 0, 2, 4],
    yticks=[-4, -2, 0, 2, 4],
)
fig.tight_layout(pad=0.5, h_pad=0.03, w_pad=0.03)
fig.savefig(figfile)
print(figfile, end='')
