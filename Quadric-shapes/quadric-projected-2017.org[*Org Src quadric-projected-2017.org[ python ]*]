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

Rcs = [0.5, 1.0, 2.0, 4.0, 8.0]
Tcs = [-2.0, -1.0, -0.5, 1e-6, 0.5, 1.0, 2.0]

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

def Rc_dash(inc, Tc, Rc):
    f = np.sqrt(1.0 + Tc*np.tan(inc)**2)
    return Rc * (1 + np.tan(inc)**2) / f / (1.0 + Rc*(f - 1.0) / Tc)

def th90(inc, Tc, Rc):
    t2i = np.tan(inc)**2
    tan_th90 = -np.sqrt(t2i*(2.0 + Tc*t2i) + (2.0 - Tc/Rc)/Rc)/t2i
    return 180.0 + np.degrees(np.arctan(tan_th90))

for Rc, lw, alpha, ls, dash in list(zip(Rcs, lws, alphas, lss, dashes))[::-1]:
    for Tc, col in list(zip(Tcs, cols))[::-1]:
        if Rc == 0.5:
            continue
        if Rc == 1.0:
            label = fr'$T_c = {Tc:.1f}$'
        else:
            label = None
        ax.plot(inc_deg, th90(inc, Tc, Rc),
                lw=1, c=col, dashes=dash, label=label)

i25, i50, i75 = [90.0 - np.degrees(np.arccos(_)) for _ in [0.25, 0.5, 0.75]]

ax.fill_betweenx([0.0, 200.0], [0.0]*2, [i25]*2, alpha=0.2, color='g')
ax.fill_betweenx([0.0, 200.0], [i50]*2, [i75]*2, alpha=0.2, color='g')
ax.legend(ncol=1, fontsize='xx-small', frameon=True)
ax.set(
    yscale='linear',
    xlim=[0.0, 90.0],
    ylim=[80.0, 160.0],
    # yticks=[1.0, 2.0, 5.0, 10.0],
    # yticklabels=['1', '2', '5', '10'],
    xlabel=r'Inclination, degrees',
    ylabel=r"Body-frame polar angle of perpendicular projected axis: $\theta_{90}$, degrees",
    xticks=[15, 30, 45, 60, 75, 90],
)        
yaxis = ax.get_yaxis()


fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
