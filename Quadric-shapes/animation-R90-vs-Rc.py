import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker
from matplotlib import animation
import seaborn as sns

fileroot = sys.argv[0].replace('.py', '')

sns.set_style('white')
sns.set_color_codes(palette='deep')
fig, ax = plt.subplots(figsize=(5, 5))

#
# Plot the background regions
#
Rc_grid = np.linspace(0.0, 10.0, 2000)
R90_T0_grid = np.sqrt(2*Rc_grid)
R90_T1_grid = np.sqrt(2*Rc_grid - 1.0)
R90_T1_grid[~np.isfinite(R90_T1_grid)] = 0.0 
ax.fill_between(Rc_grid, R90_T1_grid, R90_T0_grid, color='k', alpha=0.2)
ax.fill_between(Rc_grid, R90_T0_grid, color='k', alpha=0.1)
ax.plot(Rc_grid, R90_T0_grid, c='k', lw=0.5)
ax.axhline(1.0, lw=0.5, alpha=0.5, color='k', zorder=-1)
ax.axvline(1.0, lw=0.5, alpha=0.5, color='k', zorder=-1)
ax.plot([0.0, 10.0], [0.0, 10.0], lw=0.5, alpha=0.5, color='k', zorder=-1)

#
# Plot lines for projected tracks versus inclination
#
inc = np.linspace(0.0, 0.5*np.pi, 500, endpoint=False)
inc_deg = np.degrees(inc)

Rcs = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
Tcs = [-2.0, -1.5, -1.0, -0.5, -0.25, 1e-8, 0.25, 0.5, 1.0, 1.5, 2.0]
shapes =  ['Hyperbola']*5 + ['Parabola', 'Prolate', 'Prolate',
                             'Sphere', 'Oblate', 'Oblate', ]

n_Rc = len(Rcs)
n_Tc = len(Tcs)

lws = np.linspace(0.5, 2.0, n_Rc)
alphas = np.linspace(1.0, 0.2, n_Rc)
cols = sns.color_palette('magma', n_colors=n_Tc)
# cols = 'cgbkmry'


def Rc_prime(inc, Tc, Rc):
    f = np.sqrt(1.0 + Tc*np.tan(inc)**2)
    return Rc * (1 + np.tan(inc)**2) / f / (1.0 + Rc*(f - 1.0) / Tc)

def Tc_prime(inc, Tc):
    fsquared = 1.0 + Tc*np.tan(inc)**2
    return Tc * (1.0 + np.tan(inc)**2) / fsquared

def R90_prime(inc, Tc, Rc):
    return np.sqrt(2*Rc_prime(inc, Tc, Rc) - Tc_prime(inc, Tc))


dot_artists = {}
for Rc, lw, alpha in list(zip(Rcs, lws, alphas))[::-1]:
    for Tc, col, shape in list(zip(Tcs, cols, shapes))[::-1]:
        if Rc == 1.0:
            label = fr'{shape}: $T_c = {Tc:.1f}$'
        else:
            label = None
        ax.plot(Rc_prime(inc, Tc, Rc), R90_prime(inc, Tc, Rc),
                c=col, lw=lw, label=label)
        # Populate dict of artists that will be animated later
        dot_artists[(Rc, Tc)], = ax.plot([], [], '.', ms=10, c=col, zorder=100)

ax.legend(ncol=2, fontsize='x-small', frameon=True, loc='upper left', title='Quadric Shape')
ax.set(
    yscale='linear',
    xlim=[0.0, 4.1],
    ylim=[0.0, 4.1],
    xlabel=r"Projected planitude: $\Pi'$",
    ylabel=r"Projected alatude: $\Lambda'$",
)        
sns.despine()
fig.tight_layout()

#
# Animation of the dots
# 
def animate_dots(inclination):
    """For each quadric, update (x, y) of dot for particular `inclination`"""
    for Rc, lw, alpha in list(zip(Rcs, lws, alphas))[::-1]:
        for Tc, col, shape in list(zip(Tcs, cols, shapes))[::-1]:
            x = Rc_prime(inclination, Tc, Rc)
            y = R90_prime(inclination, Tc, Rc)
            dot_artists[(Rc, Tc)].set_data([x], [y])
    return dot_artists.values()

sini = np.linspace(0.0, 1.0, 100)
incs4anim = np.arcsin(sini)
anim = animation.FuncAnimation(fig, animate_dots, frames=incs4anim,
                               repeat_delay=100, blit=True)
moviefile = fileroot + '.mp4'
anim.save(moviefile, writer='ffmpeg', fps=30, dpi=200)
print(moviefile, end='')
