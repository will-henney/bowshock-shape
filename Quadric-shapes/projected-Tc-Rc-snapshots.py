import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker
import seaborn as sns

plotfile = sys.argv[0].replace('.py', '.png')

sns.set_style('white')
fig, axes = plt.subplots(3, 3, figsize=(9, 9), sharex=True, sharey=True)

incs_deg = 10.0*np.arange(9)

ny, nx = 65, 73
Rcs = np.linspace(0.5, 8.0, nx)
Tcs = np.linspace(-3.0, 2.0, ny)[::-1]
Rc_grid = Rcs[None, :]*np.ones_like(Tcs[:, None])
Tc_grid = Tcs[:, None]*np.ones_like(Rcs[None, :])

cols = sns.color_palette('magma', n_colors=ny)


def Rc_prime(inc, Tc, Rc):
    f = np.sqrt(1.0 + Tc*np.tan(inc)**2)
    return Rc * (1 + np.tan(inc)**2) / f / (1.0 + Rc*(f - 1.0) / Tc)

def Tc_prime(inc, Tc):
    fsquared = 1.0 + Tc*np.tan(inc)**2
    return Tc * (1.0 + np.tan(inc)**2) / fsquared

for ax, inc_deg in zip(axes.flat, incs_deg):
    inc = np.radians(inc_deg)
    Rcp = Rc_prime(inc, Tc_grid, Rc_grid).ravel()
    Tcp = Tc_prime(inc, Tc_grid).ravel()
    ax.scatter(Rcp, Tcp, c=Tc_grid.ravel(),
	       vmin=Tc_grid.min(), vmax=Tc_grid.max(),
	       edgecolors='none',
	       cmap='magma', marker='.', s=15, alpha=0.8)
    ax.axhspan(0.0, 10.0, alpha=0.1, facecolor='k', zorder=-1)
    ax.axhline(1.0, ls='--', lw=0.5, c='k', zorder=0)
    ax.axvline(1.0, ls='--', lw=0.5, c='k', zorder=0)
    ax.plot([1.0], [1.0], 'x', c='k')
    ax.text(5.5, -4.0, rf'$i = {inc_deg:.0f}^\circ$',
            bbox={'facecolor': 'w', 'alpha': 0.8, 'edgecolor': 'none'})

axes[-1, 0].set(
    yscale='linear',
    xlim=[0.0, 8.1],
    ylim=[-5.0, 2.1],
    xlabel=r"$\widetilde{R}_{c}{}'$",
    ylabel=r"$T_c{}'$",
)        

fig.tight_layout()
fig.savefig(plotfile, dpi=300)
print(plotfile, end='')
