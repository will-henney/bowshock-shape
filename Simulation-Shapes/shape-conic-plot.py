import sys
import json
import numpy as np
from astropy.coordinates import Longitude
from matplotlib import pyplot as plt
import seaborn as sns
sys.path.append("../Dust-wave")
import bow_projection as bp

try:
    prefix = sys.argv[1]
except:
    print(f"Usage: {sys.argv[0]} PREFIX")

plotfile = sys.argv[0].replace('.py', f'-{prefix}.pdf')

def compensate(R, theta):
    """Compensated inversion of R(theta)"""
    return 1.0/R - 0.5*(1 + np.cos(theta))

def load_R_th(arc_prefix):
    jfile = f'{arc_prefix}-arcdata.json'
    data = json.load(open(jfile))
    R0 = np.array(data['outer']['R0'])
    R = np.array(data['outer']['R'])
    th = Longitude(data['outer']['theta'], unit='deg')
    th += Longitude(data['outer']['PA0'], unit='deg')
    return th.rad, R/R0


sns.set_style('ticks')
fig, ax = plt.subplots()

# Plot confocal parabola
ax.axhline(0.0, ls='-', c='k', lw=0.5)

# Plot wilkinoid
mugrid = np.linspace(-1.0, 1.0, 200)
thgrid = np.arccos(mugrid)
ax.plot(mugrid, compensate(bp.wilkinoid_R_theta(thgrid), thgrid),
        '-', c='k', lw=1.5)

# Plot cantoids
for beta in 0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 0.0003, 0.0001:
    ax.plot(mugrid, compensate(bp.cantoid_R_theta(thgrid, beta), thgrid),
            '-', c='k', lw=0.5)

# Fill in forbidden zone
ax.fill_between(mugrid, -0.5*(1.0 + mugrid), -1.0, color='k', alpha=0.4)

# Plot traced arcs
try:
    th, R = load_R_th(prefix + '-CD')
    ax.plot(np.cos(th), compensate(R, th), '.', alpha=0.6, label='CD')
    th, R = load_R_th(prefix + '-BS')
    ax.plot(np.cos(th), compensate(R, th), '.', alpha=0.6, label='BS')
except:
    th, R = load_R_th(prefix)
    ax.plot(np.cos(th), compensate(R, th), '.', alpha=0.6, label=prefix)



ax.legend(title=prefix)

ax.set(
    xlim=[-1.02, 1.02],
    ylim=[-0.155, 0.155],
    xlabel=r"$\cos \,\theta$",
    ylabel=r"$(R_{0} / R) - 0.5 (1 + \cos \,\theta) $",
)
sns.despine(trim=True)
fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
