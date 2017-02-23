import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker
import seaborn as sns

plotfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('white')
fig, ax = plt.subplots(figsize=(5, 5))

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

inc = np.linspace(0.0, 0.5*np.pi, 500, endpoint=False)
inc_deg = np.degrees(inc)

Rcs = [0.5, 1.0, 2.0, 4.0, 8.0]
Tcs = [-2.0, -1.0, -0.5, 1e-8, 0.5, 1.0, 2.0]

n_Rc = len(Rcs)
n_Tc = len(Tcs)

lws = np.linspace(1.0, 2.0, n_Rc)
dash_solid = []
dash_dashed = [3, 2]
dash_dotted = [1, 2]
dash_dot_dashed = [1, 2, 4, 2]
dash_triple_dot_dashed = [1, 2, 1, 2, 1, 2, 4, 2]
dashes = [dash_triple_dot_dashed, dash_solid,
          dash_dashed, dash_dotted, dash_dot_dashed]

lss = ['-.', '-', '--', ':', '-.']
alphas = np.linspace(1.0, 0.2, n_Rc)
cols = sns.color_palette('magma', n_colors=n_Tc)


def Rc_prime(inc, Tc, Rc):
    f = np.sqrt(1.0 + Tc*np.tan(inc)**2)
    return Rc * (1 + np.tan(inc)**2) / f / (1.0 + Rc*(f - 1.0) / Tc)

def Tc_prime(inc, Tc):
    fsquared = 1.0 + Tc*np.tan(inc)**2
    return Tc * (1.0 + np.tan(inc)**2) / fsquared

def R90_prime(inc, Tc, Rc):
    return np.sqrt(2*Rc_prime(inc, Tc, Rc) - Tc_prime(inc, Tc))


for Rc, lw, alpha, dash in list(zip(Rcs, lws, alphas, dashes))[::-1]:
    for Tc, col in list(zip(Tcs, cols))[::-1]:
        if Rc == 1.0:
            label = fr'$T_c = {Tc:.1f}$'
        else:
            label = None
        ax.plot(Rc_prime(inc, Tc, Rc), R90_prime(inc, Tc, Rc),
                c=col, dashes=dash, label=label)
        # ax.plot(Rc_dash(inc, Tc, Rc), Tc_dash(inc, Tc), '.', alpha=0.1, ms=4,
        #         c=col, label=label)
        ax.plot([Rc_prime(0, Tc, Rc)], [R90_prime(0, Tc, Rc)], 'o', c=col)

ax.legend(ncol=1, fontsize='xx-small', frameon=True)
ax.set(
    yscale='linear',
    xlim=[0.0, 8.1],
    ylim=[0.0, 8.1],
    xlabel=r"Projected dimensionless radius of curvature: $\widetilde{R}_{c}{}'$",
    ylabel=r"Projected dimensionless perpendicular radius: $\widetilde{R}_{90}{}'$",
)        

fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
