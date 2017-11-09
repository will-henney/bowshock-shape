import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker
import seaborn as sns
import conic_parameters

plotfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('ticks')
sns.set_color_codes(palette='deep')
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

# Plot our curves
betas = np.logspace(-5.0, np.log10(0.5), 100)
ks = [0.0, 0.5, 3.0, 8.0]
labels = ["Cantoid"] + [f"Ancantoid, $k = {k:.1f}$" for k in ks[1:]]
bbetas = [1e-8, 0.01, 0.1]
mss = [6, 4, 3]
colors = 'krmb'
ax.plot(5./3., np.sqrt(3.0), '+', color='white', ms=12)
for k, color, label in zip(ks, colors, labels):
    xi = 2.0 / (k + 2.0)
    Pi = conic_parameters.A(betas, xi)
    Lambda = conic_parameters.B(betas, xi)
    ax.plot(Pi, Lambda, label=label, color=color, alpha=0.7)
    for bb, ms in zip(bbetas, mss):
        ax.plot(conic_parameters.A(bb, xi), conic_parameters.B(bb, xi), 'o',
                label="_nolabel_", color=color, ms=ms)

    # Test accuracy of approximate solution at beta = 0.3
    ax.plot(conic_parameters.A(0.5, xi),
            conic_parameters.B(0.5, xi,
                               th1_90=conic_parameters.th1_90_exact),
            's',
            label="_nolabel_", color=color, ms=2)
ax.legend(ncol=1, fontsize='small', frameon=True, loc='upper left', title='Shape')
ax.set(
    yscale='linear',
    xlim=[0.0, 5.75],
    ylim=[0.0, 5.75],
    xlabel=r"True planitude: $\Pi$",
    ylabel=r"True alatude: $\Lambda$",
)
ax.set_aspect('equal')
sns.despine()
fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
