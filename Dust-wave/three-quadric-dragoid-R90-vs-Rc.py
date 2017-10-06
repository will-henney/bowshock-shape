import sys
import json
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker
import seaborn as sns

def Rc_prime(inc, Tc, Rc):
    f = np.sqrt(1.0 + Tc*np.tan(inc)**2)
    return Rc * (1 + np.tan(inc)**2) / f / (1.0 + Rc*(f - 1.0) / Tc)

def Tc_prime(inc, Tc):
    fsquared = 1.0 + Tc*np.tan(inc)**2
    return Tc * (1.0 + np.tan(inc)**2) / fsquared

def R90_prime(inc, Tc, Rc):
    return np.sqrt(2*Rc_prime(inc, Tc, Rc) - Tc_prime(inc, Tc))

plotfile = sys.argv[0].replace('.py', '.pdf')

alldata = json.load(open('dust-wave-fitdata.json'))


sns.set_style('white')
sns.set_color_codes('dark')

fig, ax = plt.subplots(figsize=(5, 5))


left_annotate_pars = dict(xytext=(-5, 5), ha='right', va='bottom')
right_annotate_pars = dict(xytext=(5, -5), ha='left', va='top')


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

colors = 'bmgr'

for (alpha, data), color in zip(alldata.items(), colors):
    # Parameters for head conic
    R0_h = 1.0
    T_h = data['head']['T']
    tilde_Rc_h = data['head']['Rc']
    R90_h = data['head']['R90']
    ax.plot([tilde_Rc_h], [R90_h], 'o', color=color)
    ax.plot(Rc_prime(inc, T_h, tilde_Rc_h),
            R90_prime(inc, T_h, tilde_Rc_h), '--', color=color)

# Put a cross at the Wilkinoid coordinates: [5/3, sqrt(3)]
ax.plot([5./3.], [np.sqrt(3.0)], '+', c='w', ms=10, alpha=1.0)

ax.legend(ncol=1, fontsize='small', frameon=True)
ax.set(
    yscale='linear',
    xscale='linear',
    xlim=[0.0, 3.1],
    ylim=[0.0, 3.1],
    xlabel=r"Projected dimensionless radius of curvature: $\widetilde{R}_{c}{}'$",
    ylabel=r"Projected dimensionless perpendicular radius: $\widetilde{R}_{90}{}'$",
)        

sns.despine()
fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
