import sys
import scanf
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import seaborn as sns
import bow_projection as bp
import standing_wave

try:
    BASE_SHAPE_ID = sys.argv[1]
except:
    sys.exit(f"Usage: {sys.argv[0]} BASE_SHAPE_ID")


fileroot = sys.argv[0].replace('.py', f'-{BASE_SHAPE_ID}')

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
fig, axes = plt.subplots(1, 3, figsize=(8, 4))

# Different (amplitude, wavenumber) are different panels
amplitudes = [0.1, 0.05, 0.02]
wavenumbers = [1.0, 2.0, 5.0]

# All are for zero inclination
inclination = 0.0

th_inf = bp.theta_infinity(base_shape)

line_artists = {}
for amplitude, wavenumber, ax in zip(amplitudes, wavenumbers, axes.flat):
    wave_label = fr"$A = {amplitude}$, $N = {wavenumber}$"
    # artist that will get animated for the different phases
    line_artists[(amplitude, wavenumber)], = ax.plot([], [],
                                                     lw=2, alpha=0.7,
                                                     color='r', label='perturbed')
    # plot the bsse shape
    th = np.linspace(0.0, th_inf, 301)
    xp, yp = bp.xyprime_t(th, inclination, base_shape)
    m = np.isfinite(xp) & np.isfinite(yp)
    xxp = np.concatenate((xp[m][::-1], xp[m]))
    yyp = np.concatenate((-yp[m][::-1], yp[m]))
    ax.plot(xxp, yyp, 'k', lw=1, label='base')
    # 
    ax.plot([0], [0], 'o', color='k')
    ax.axhline(0, color='k', lw=0.8)
    ax.legend(title=wave_label, ncol=1, loc="upper right")
    ax.set(
        xlabel=r"$x / R_0$",
        ylabel=r"$y / R_0$",
        xlim=[-2.5, 1.5],
        ylim=[-0.5, 4.8],
        xticks=[-2, -1, 0, 1],
        yticks=[0, 1, 2, 3],
    )
    ax.set_aspect('equal', adjustable='box')


def animate_bows(phase):
    for amplitude, wavenumber, ax in zip(amplitudes, wavenumbers, axes.flat):
        shape = standing_wave.StandingWave(base_shape,
                                           amplitude=amplitude,
                                           wavenumber=wavenumber)
        th_inf = bp.theta_infinity(shape)
        shape.phase = phase
        th0, th90 = bp.theta_0_90(inclination, shape)
        th = np.linspace(th0, th_inf, 301)
        xp, yp = bp.xyprime_t(th, inclination, shape)
        m = np.isfinite(xp) & np.isfinite(yp)
        if m.sum() == 0:
            # Case of no tangent line at all at this inclination
            continue
        xxp = np.concatenate((xp[m][::-1], xp[m]))
        yyp = np.concatenate((-yp[m][::-1], yp[m]))
        radii = bp.characteristic_radii_projected(inclination, shape)
        R0p = radii['R_0 prime']
        R0p = 1.0
        x = xxp/R0p
        y = yyp/R0p
        line_artists[(amplitude, wavenumber)].set_data(x, y)
    return line_artists.values()

sns.despine(bottom=True)
fig.tight_layout()

anim = animation.FuncAnimation(fig, animate_bows,
                               frames=np.linspace(0, 1.0, 50),
                               blit=True)
moviefile = fileroot + '.mp4'
anim.save(moviefile, writer='ffmpeg', fps=60, dpi=200)
print(moviefile, end='')
