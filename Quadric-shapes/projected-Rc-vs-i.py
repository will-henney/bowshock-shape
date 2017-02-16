import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker
import seaborn as sns

plotfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('white')
fig, ax = plt.subplots(figsize=(5, 5))

inc = np.linspace(0.0, 0.5*np.pi, 500)
inc_deg = np.degrees(inc)

Rcs = [1.0, 2.0, 4.0, 8.0]
Tcs = [-2.0, -1.0, -0.5, 1e-3, 0.5, 1.0, 2.0]

n_Rc = len(Rcs)
n_Tc = len(Tcs)

lws = np.linspace(1.0, 2.0, n_Rc)
lss = ['-', '--', ':', '-.']
alphas = np.linspace(1.0, 0.2, n_Rc)
cols = sns.color_palette('magma', n_colors=n_Tc)


def Rc_dash(inc, Tc, Rc):
    f = np.sqrt(1.0 + Tc*np.tan(inc)**2)
    return Rc * (1 + np.tan(inc)**2) / f / (1.0 + Rc*(f - 1.0) / Tc)

for Rc, lw, alpha, ls in list(zip(Rcs, lws, alphas, lss))[::-1]:
    for Tc, col in list(zip(Tcs, cols))[::-1]:
        if Rc == 1.0:
            label = fr'$T_c = {Tc:.1f}$'
        else:
            label = None
        ax.plot(inc_deg, Rc_dash(inc, Tc, Rc),
                c=col, ls=ls, label=label)

i25, i50, i75 = [90.0 - np.degrees(np.arccos(_)) for _ in [0.25, 0.5, 0.75]]

ax.fill_betweenx([0.0, 100.0], [0.0]*2, [i25]*2, alpha=0.2, color='g')
ax.fill_betweenx([0.0, 100.0], [i50]*2, [i75]*2, alpha=0.2, color='g')
ax.legend(ncol=1, fontsize='xx-small', frameon=True)
ax.set(
    yscale='linear',
    xlim=[0.0, 90.0],
    ylim=[0.0, 10.0],
    # yticks=[1.0, 2.0, 5.0, 10.0],
    # yticklabels=['1', '2', '5', '10'],
    xlabel=r'Inclination, degrees',
    ylabel=r"Projected dimensionless radius of curvature: $\widetilde{R}_{c}{}'$",
    xticks=[15, 30, 45, 60, 75, 90],
)        
yaxis = ax.get_yaxis()

# yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=2.0))
# yaxis.set_major_formatter(matplotlib.ticker.LogFormatter())

yaxis.set_major_locator(matplotlib.ticker.FixedLocator([1.0, 2.0, 4.0, 8.0]))

fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
