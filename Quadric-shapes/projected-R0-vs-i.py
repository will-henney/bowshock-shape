import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

plotfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('white')
fig, ax = plt.subplots(figsize=(5, 5))

inc = np.linspace(0.0, 0.5*np.pi, 5000)
inc_deg = np.degrees(inc)

Rcs = [0.5, 1.0, 2.0, 4.0, 8.0]
Tcs = [-2.0, -1.0, -0.5, 1e-3, 0.5, 1.0, 2.0]

n_Rc = len(Rcs)
n_Tc = len(Tcs)

lws = np.linspace(1.0, 2.0, n_Rc)
alphas = np.linspace(1.0, 0.2, n_Rc)
cols = sns.color_palette('magma', n_colors=n_Tc)
lss = ['-.', '-', '--', ':', '-.']

dash_solid = []
dash_dashed = [3, 2]
dash_dotted = [1, 2]
dash_dot_dashed = [1, 2, 4, 2]
dash_triple_dot_dashed = [1, 2, 1, 2, 1, 2, 4, 2]
dashes = [dash_triple_dot_dashed, dash_solid,
          dash_dashed, dash_dotted, dash_dot_dashed]

def qratio(inc, Tc, Rc):
    f = np.sqrt(1.0 + Tc*np.tan(inc)**2)
    return (1.0 + Rc*(f - 1.0) / Tc)*np.cos(inc)

for Rc, lw, alpha, ls, dash in list(zip(Rcs, lws, alphas, lss, dashes))[::-1]:
    for Tc, col in list(zip(Tcs, cols))[::-1]:
        if Rc == 1.0:
            label = fr'$T_c = {Tc:.1f}$'
        else:
            label = None
        ax.plot(inc_deg, qratio(inc, Tc, Rc),
                c=col, dashes=dash, label=label)

i25, i50, i75 = [90.0 - np.degrees(np.arccos(_)) for _ in [0.25, 0.5, 0.75]]

ax.fill_betweenx([0.0, 100.0], [0.0]*2, [i25]*2, alpha=0.2, color='g')
ax.fill_betweenx([0.0, 100.0], [i50]*2, [i75]*2, alpha=0.2, color='g')

ax.plot(inc_deg, np.cos(inc), c='b')
ax.annotate("Variation in projected\nseparation " + r"$D'/D$",
            xy=(50.0, np.cos(np.radians(50.0))), xycoords='data',
            xytext=(-20, -10), textcoords='offset points',
            ha='right', va='top',
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"),
            fontsize='x-small',)

ax.legend(ncol=1, fontsize='xx-small', frameon=True)
ax.set(
    yscale='linear',
    xlim=[0.0, 90.0],
    ylim=[0.0, 5.5],
    xlabel=r'Inclination, degrees',
    ylabel=r"Variation in projected stand-off distance: $R_{0}' / R_{0}$",
    xticks=[15, 30, 45, 60, 75, 90],
)        
fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
