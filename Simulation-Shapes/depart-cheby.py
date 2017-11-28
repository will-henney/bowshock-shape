import sys
import json
import numpy as np
from astropy.modeling import models, fitting
from astropy.coordinates import Longitude
from matplotlib import pyplot as plt
import seaborn as sns
sys.path.append("../Dust-wave")
import bow_projection as bp

try:
    prefix = sys.argv[1]
except:
    print(f"Usage: {sys.argv[0]} PREFIX [CHEBY_DEGREE] [EXTRAP_DEGREE]")

try:
    cheby_degree = int(sys.argv[2])
except:
    cheby_degree = 10

try:
    extrap_degree = int(sys.argv[3])
except:
    extrap_degree = 2

plotfile = sys.argv[0].replace('.py', f'-{prefix}.pdf')

def load_R_th(arc_prefix):
    jfile = f'{arc_prefix}-arcdata.json'
    data = json.load(open(jfile))
    R0 = np.array(data['outer']['R0'])
    R = np.array(data['outer']['R'])
    th = Longitude(data['outer']['theta'], unit='deg')
    th += Longitude(data['outer']['PA0'], unit='deg')
    return th.rad, R/R0


def departure(R, theta):
    """Parabolic departure of R(theta)"""
    return 1.0/R - 0.5*(1 + np.cos(theta))

def extrapolate(mu, Delta, mu0=-0.5, force_open=False, deg=2):
    def factor(mu):
        if force_open:
            return np.abs(-1.0 - mu)**0.5
        else:
            return 1.0

    # Only fit mu < mu0
    mask = mu <= mu0
    p = np.poly1d(np.polyfit(mu[mask], Delta[mask]/factor(mu[mask]), deg=deg))
    mu_x = np.linspace(-1.0, mu0)
    return mu_x, factor(mu_x)*p(mu_x)


sns.set_style('ticks')
sns.set_color_codes('deep')
fig, ax = plt.subplots(figsize=(4, 4))

# Plot x=0, y=0 axes
ax.axhline(0.0, ls=':', c='k', lw=0.5)
ax.axvline(0.0, ls=':', c='k', lw=0.5)

# Plot wilkinoid
mugrid = np.linspace(-1.0, 1.0, 200)
thgrid = np.arccos(mugrid)
ax.plot(mugrid, departure(bp.wilkinoid_R_theta(thgrid), thgrid),
        '-', c='k', lw=1.5)

# Plot cantoids
for beta in 0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 0.0003, 0.0001:
    ax.plot(mugrid, departure(bp.cantoid_R_theta(thgrid, beta), thgrid),
            '-', c='k', lw=0.5)

# Fill in forbidden zone
# ax.fill_between(mugrid, -0.5*(1.0 + mugrid), -1.0, color='k', alpha=0.4)

# Plot traced arcs
th, R = load_R_th(prefix)
Delta = departure(R, th)
mu = np.cos(th)
T = models.Chebyshev1D(degree=cheby_degree)
fitter = fitting.LevMarLSQFitter()
T = fitter(T, mu, Delta)
ax.plot(mu, T(mu), '-', alpha=0.5, color='r', lw=4, label='_nolabel_')

mux, Deltax = extrapolate(mu, Delta, deg=extrap_degree, force_open=False)
ax.plot(mux, Deltax, '--', color='r', label='_nolabel_')
mux, Deltax = extrapolate(mu, Delta, deg=extrap_degree, force_open=True)
ax.plot(mux, Deltax, '-', color='r',  label='_nolabel_')

ax.plot(mu, Delta, '.', color='b', alpha=0.8, label=prefix)

title = "MHD simulation" if "MHD" in prefix else "HD simulation"
ax.text(0.5, 0.1, title, ha='center', va='bottom')

ax.set(
    xlim=[-1.05, 1.05],
    ylim=[-0.155, 0.155],
    xlabel=r"$\cos \,\theta$",
    ylabel=r"Parabolic departure function, $\Delta(\cos\theta)$",
)
sns.despine(trim=True)
fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
