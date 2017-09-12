import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import FixedLocator
import seaborn as sns

plotfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('white')
sns.set_color_codes(palette='deep')
fig, ax = plt.subplots(figsize=(5, 5))

inc = np.linspace(0.0, 0.5*np.pi, 500, endpoint=False)
inc_deg = np.degrees(inc)

Rcs = [0.5, 1.0, 2.0, 4.0, 8.0]
Tcs = [-2.0, -1.0, -0.5, 1e-12, 0.5, 1.0, 2.0]
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

def Rc_prime(inc, Tc, Rc):
    f = np.sqrt(1.0 + Tc*np.tan(inc)**2)
    return Rc * (1 + np.tan(inc)**2) / f / (1.0 + Rc*(f - 1.0) / Tc)

def Tc_prime(inc, Tc):
    fsquared = 1.0 + Tc*np.tan(inc)**2
    return Tc * (1.0 + np.tan(inc)**2) / fsquared

def R90_prime(inc, Tc, Rc):
    return np.sqrt(2*Rc_prime(inc, Tc, Rc) - Tc_prime(inc, Tc))

for Rc, lw, alpha, ls, dash in list(zip(Rcs, lws, alphas, lss, dashes))[::-1]:
    for Tc, col, shape in list(zip(Tcs, cols, shapes))[::-1]:
        if Rc == 1.0:
            label = fr'{shape}: $T_c = {Tc:.1f}$'
        else:
            label = None
        ax.plot(inc_deg, R90_prime(inc, Tc, Rc),
                c=col, lw=lw, alpha=alpha, label=label)

i25, i50, i75 = [90.0 - np.degrees(np.arccos(_)) for _ in [0.25, 0.5, 0.75]]

ax.fill_betweenx([0.0, 100.0], [0.0]*2, [i25]*2, alpha=0.1, color='k')
ax.fill_betweenx([0.0, 100.0], [i50]*2, [i75]*2, alpha=0.1, color='k')
ax.legend(ncol=1, fontsize='small', frameon=True, borderaxespad=0, title='Quadric Shape')
ax.set(
    yscale='linear',
    xlim=[0.0, 90.0],
    ylim=[0.0, 5.0],
    # yticks=[1.0, 2.0, 5.0, 10.0],
    # yticklabels=['1', '2', '5', '10'],
    xlabel=r'Inclination, degrees',
    ylabel=r"Projected dimensionless perpendicular radius: $\widetilde{R}_{90}{}'$",
    xticks=[15, 30, 45, 60, 75, 90],
)        
yaxis = ax.get_yaxis()

# yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=2.0))
# yaxis.set_major_formatter(matplotlib.ticker.LogFormatter())

yaxis.set_major_locator(FixedLocator([1.0, 2.0, 3.0, 4.0]))
sns.despine()

fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
