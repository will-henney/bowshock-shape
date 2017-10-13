import sys
import scanf
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import bow_projection as bp
import ancantoid_shape
import dragoid_shape
import standing_wave

try:
    AMPLITUDE = float(sys.argv[1])
    WAVENUMBER = float(sys.argv[2])
    BASE_SHAPE_ID = sys.argv[3]
except:
    sys.exit(f"Usage: {sys.argv[0]} AMPLITUDE WAVENUMBER BASE_SHAPE_ID")


figfile = sys.argv[0].replace(
    '.py', f'-A{int(100*AMPLITUDE):03d}-N{int(10*WAVENUMBER):02d}'
    f'-{BASE_SHAPE_ID}.pdf')

# Choose which base shape according to command-line argument and parse
# out the shape parameters if any
if BASE_SHAPE_ID == "paraboloid":
    base_shape = bp.paraboloid_R_theta
    shape_label = "Paraboloid"
elif BASE_SHAPE_ID == "wilkinoid":
    base_shape = bp.wilkinoid_R_theta
    shape_label = "Wilkinoid"
elif BASE_SHAPE_ID.startswith("cantoid"):
    ibeta, = scanf.scanf("cantoid-beta%d", BASE_SHAPE_ID)
    beta = ibeta / 100000
    base_shape = bp.Spline_R_theta_from_function(
        ngrid=1000, shape_func=bp.cantoid_R_theta, shape_func_pars=(beta,))
    shape_label = rf"Cantoid $\beta = {beta}$"
elif BASE_SHAPE_ID.startswith("ancantoid"):
    ixi, ibeta = scanf.scanf("ancantoid-xi%d-beta%d", BASE_SHAPE_ID)
    xi, beta = ixi / 100, ibeta / 100000
    base_shape = ancantoid_shape.Ancantoid(xi=xi, beta=beta, n=301)
    shape_label = rf"Ancantoid $\xi = {xi:.1f}$, $\beta = {beta}$"
elif BASE_SHAPE_ID.startswith("dragoid"):
    ialpha, = scanf.scanf("dragoid-alpha%d", BASE_SHAPE_ID)
    alpha = ialpha / 100
    base_shape = dragoid_shape.Dragoid(alpha=alpha)
    shape_label = rf"Dragoid $\alpha_\mathrm{{drag}} = {alpha:.2f}$"


sns.set_style('ticks')
fig, axes = plt.subplots(1, 3, figsize=(9, 4))

inclinations = [0, 15, 30, 45, 60, 75]
linewidths = [2.4, 2.0, 1.6, 1.2, 0.8, 0.4]
colors = sns.color_palette(palette="magma_r", n_colors=len(inclinations))


shape = standing_wave.StandingWave(
    base_shape, amplitude=AMPLITUDE, wavenumber=WAVENUMBER)

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
        xlim=[-2.5, 1.2],
        ylim=[-3, 3],
    )
    ax.set_aspect('equal', adjustable='box')

sns.despine()
fig.tight_layout()
fig.savefig(figfile)
print(figfile, end='')
