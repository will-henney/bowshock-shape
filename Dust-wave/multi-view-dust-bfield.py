import sys
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import seaborn as sns

try:
    ANGLE = int(sys.argv[1])
except:
    sys.exit(f"Usage: {sys.argv[0]} ANGLE")


figfile = sys.argv[0].replace('.py', f'-{ANGLE:02d}.pdf')
fig, axes = plt.subplots(3, 4, sharex=True, sharey=True, figsize=(8, 6.1))

xmin, xmax, ymin, ymax = -5.0, 5.0, -5.0, 5.0
for j, inc in enumerate([7, 22, 39, 61]):
    for i, azi in enumerate([15, 45, 75]):
        fitsfile = f"dust-bfield-b{ANGLE:02d}-cube-F-{inc:03d}-{azi:03d}.fits"
        hdu, = fits.open(fitsfile)

        ax = axes[i, j]
        ax.imshow(np.sqrt(hdu.data),
                  origin='lower', extent=[xmin, xmax, ymin, ymax])
        thBtext = fr"$\theta_B = {ANGLE}^\circ$"
        inctext = fr"$i = {inc}^\circ$"
        azitext = fr"$\phi = {azi}^\circ$"
        ax.text(-4.2, -4.2, inctext, color="w")
        ax.text(4.2, -4.2, azitext, color="w", ha='right')
        ax.text(4.2, 4.2, thBtext, color="w", ha='right', va='top')
        ax.plot([0], [0], '+', color='orange')
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
