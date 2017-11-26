import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

plotfile = sys.argv[0].replace('.py', '.pdf')

def departure(R, theta):
    """Departure function of R(theta) from parabola"""
    return 1.0/R - 0.5*(1 + np.cos(theta))


def xy_conic_t(t, Pi, Lambda):
    """Parametric equation for ellipse/hyperbola"""
    Q = 2*Pi - Lambda**2
    a = Pi / np.abs(Q)
    b = a*np.sqrt(np.abs(Q))
    sigma = np.sign(Q)
    x0 = 1.0 - sigma*a
    if sigma > 0:
        C, S = np.cos, np.sin
    else:
        C, S = np.cosh, np.sinh
    x = x0 + sigma*a*C(t)
    y = b*S(t)
    return x, y

def R_th_from_xy(x, y):
    R = np.hypot(x, y)
    th = np.arctan2(y, x)
    return R, th

sns.set_style('ticks')
sns.set_color_codes('deep')
fig, ax = plt.subplots(figsize=(4, 4))

# Show theta = 90 
ax.axvline(0.0, ls=':', c='k', lw=0.5)

t = np.linspace(0.0, np.pi, 500)
colors = 'rky'
lws = 0.5, 1.0, 1.5, 2.0
Tudes = [1.5, 1.9999, 2.66]


for Pi, color in zip(Tudes, colors):
    for Lambda, lw in zip(Tudes, lws):
        x, y = xy_conic_t(t, Pi, Lambda)
        R, th =  R_th_from_xy(x, y)
        Delta = departure(R, th)
        label = fr'$\Pi = {Pi:.1f}$, $\Lambda = {Lambda:.1f}$'
        ax.plot(np.cos(th), Delta, c=color, lw=lw, label=label)

# Fill in forbidden zone
# mugrid = np.linspace(-1.0, 1.0)
# ax.fill_between(mugrid, -0.5*(1.0 + mugrid), -1.0, color='k', alpha=0.4)

ax.annotate(r"$\Pi = \frac{3}{2}$", (0.8, 0.05), color='r', va='center')
ax.annotate(r"$\Pi = 2$", (1.01, 0.0), color='k', va='center')
ax.annotate(r"$\Pi = \frac{8}{3}$", (0.8, -0.05), color='y', va='center')

ax.annotate(r"$\Lambda = \frac{3}{2}$", (0.05, 0.18), va='bottom')
ax.annotate(r"$\Lambda = 2$", (0.0, 0.02), va='bottom', ha='center')
ax.annotate(r"$\Lambda = \frac{8}{3}$", (0.05, -0.14), va='top')

ax.annotate(r"Apex", (1.0, 0.1), fontsize='small', va='bottom', ha='center')
ax.annotate(r"Wing", (0.0, 0.3), fontsize='small', va='bottom', ha='center')
ax.annotate(r"Far wing", (-0.6, 0.05), fontsize='small', va='bottom', ha='center')

ax.set(
    xlim=[-1.02, 1.14],
    ylim=[-0.255, 0.455],
    xlabel=r"$\cos \,\theta$",
    ylabel=r"Departure function, $\Delta(\cos\theta)$",
)
sns.despine(trim=True)
fig.tight_layout()
fig.savefig(plotfile)

print(plotfile, end='')
