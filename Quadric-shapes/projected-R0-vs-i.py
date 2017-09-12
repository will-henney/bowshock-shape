import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns

plotfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('white')
sns.set_color_codes(palette='deep')
fig, ax = plt.subplots(figsize=(5, 5))

inc = np.linspace(0.0, 0.5*np.pi, 5000)
inc_deg = np.degrees(inc)

Rcs = [0.5, 1.0, 2.0, 4.0, 8.0]
thlabs = [44, 50, 42, 33, 22]
Tcs = [-2.0, -1.0, -0.5, 1e-3, 0.5, 1.0, 2.0]
shapes =  ['Hyperbola']*3 + ['Parabola', 'Prolate', 'Sphere', 'Oblate', ]

n_Rc = len(Rcs)
n_Tc = len(Tcs)

lws = np.linspace(0.5, 2.0, n_Rc)
alphas = np.linspace(1.0, 1.0, n_Rc)**1.5
cols = sns.color_palette('magma', n_colors=n_Tc)
cols = 'cgbkmry'
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

for Rc, lw, alpha, ls, dash, thlab in list(zip(Rcs, lws, alphas, lss, dashes, thlabs))[::-1]:
    for Tc, col, shape in list(zip(Tcs, cols, shapes))[::-1]:
        if Rc == 1.0:
            label = fr'{shape}: $T_c = {Tc:.1f}$'
        else:
            label = None
        ax.plot(inc_deg, qratio(inc, Tc, Rc),
                c=col, lw=lw, alpha=alpha, label=label)
        # ax.plot(inc_deg, qratio(inc, Tc, Rc),
        #         c=col, dashes=dash, label=label)
    ax.text(thlab, qratio(np.radians(thlab), 0.5, Rc),
            r'$\widetilde{R}_{c}{} = ' + f'{Rc:.1f}$',
            ha='center', va='center', fontsize='x-small',
            bbox={'facecolor': 'white', 'alpha': 0.7, 'pad': 0.05, 'ec': 'none'})

i25, i50, i75 = [90.0 - np.degrees(np.arccos(_)) for _ in [0.25, 0.5, 0.75]]

ax.fill_betweenx([0.0, 100.0], [0.0]*2, [i25]*2, alpha=0.1, color='k')
ax.fill_betweenx([0.0, 100.0], [i50]*2, [i75]*2, alpha=0.1, color='k')

ax.plot(inc_deg, np.cos(inc), c='b')
ax.annotate("Variation in projected\nseparation " + r"$D'/D$",
            xy=(50.0, np.cos(np.radians(50.0))), xycoords='data',
            xytext=(-20, -10), textcoords='offset points',
            ha='right', va='top',
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"),
            fontsize='small',)

ax.legend(ncol=1, fontsize='small', frameon=True, title='Quadric Shape')
ax.set(
    yscale='log',
    xlim=[0.0, 90.0],
    # ylim=[0.0, 5.5],
    ylim=[0.05, 50],
    xlabel=r'Inclination, degrees',
    ylabel=r"Variation in projected stand-off distance: $R_{0}' / R_{0}$",
    xticks=[15, 30, 45, 60, 75, 90],
)
ax.yaxis.set_major_formatter( FormatStrFormatter('%.1f') )

sns.despine(trim=False)
fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
