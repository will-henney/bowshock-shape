import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import bow_projection as bp
import standing_wave

try:
    amplitude = float(sys.argv[1])
    wavenumber = float(sys.argv[2])
except:
    sys.exit(f"Usage: {sys.argv[0]} AMPLITUDE WAVENUMBER")


figfile = sys.argv[0].replace(
    '.py', f'-A{int(100*amplitude):03d}-N{int(10*wavenumber):02d}.pdf')

sns.set_style('ticks')
fig, axes = plt.subplots(2, 2, figsize=(8, 8))

inclinations = [0, 15, 30, 45, 60, 75]
linewidths = [2.4, 2.0, 1.6, 1.2, 0.8, 0.4]
colors = sns.color_palette(palette="magma_r", n_colors=len(inclinations))


shape = standing_wave.StandingWave(
    bp.wilkinoid_R_theta, amplitude=amplitude, wavenumber=wavenumber)

phases = np.linspace(0.0, 0.5, len(axes.flat))
for phase, ax in zip(phases, axes.flat):
    shape.phase = phase
    label = f"Phase = ${phase:.2f}$"
    th_inf = bp.theta_infinity(shape)
    for inc_dg, color, lw in zip(inclinations, colors, linewidths):
        inc = np.radians(inc_dg)
        th0, th90 = bp.theta_0_90(inc, shape)
        th = np.linspace(th0, th_inf, 301)
        xp, yp = bp.xyprime_t(th, inc, shape)
        m = np.isfinite(xp) & np.isfinite(yp)
        if m.sum() == 0:
            # Case of no tangent line at all at this inclination
            continue
        xxp = np.concatenate((xp[m][::-1], xp[m]))
        yyp = np.concatenate((-yp[m][::-1], yp[m]))
        radii = bp.characteristic_radii_projected(inc, shape)
        R0p = radii['R_0 prime']
        ax.plot(xxp/R0p, yyp/R0p,
                label=fr"$i = {inc_dg:d}^\circ$",
                color=color, lw=1.5*lw)

    ax.plot([0], [0], 'o', color='k')

    ax.legend(title=label, ncol=1, loc="center left")
    ax.set(
        xlabel=r"$x' / R_0'$",
        ylabel=r"$y' / R_0'$",
        xlim=[-7, 3],
        ylim=[-5, 5],
    )
    ax.set_aspect('equal', adjustable='box')

sns.despine()
fig.tight_layout()
fig.savefig(figfile)
print(figfile, end='')
