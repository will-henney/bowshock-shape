import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

figfile = sys.argv[0].replace('.py', '.pdf')
NTH = 4001
sns.set_style('ticks')
sns.set_color_codes('dark')
fig, ax = plt.subplots(figsize=(4, 4))
blist = np.linspace(0.0, 6.0) + 0.01
thmlist = np.arccos(1./np.sqrt(1.0 + 4.0*blist**2))
for thm in thmlist:
    epsilon = 1./np.cos(thm)
    theta = np.linspace(0.001, min(np.pi, 2*thm - 0.001), NTH)
    im = np.argmin(np.abs(theta - thm))
    r = 0.5*(epsilon**2 - 1)/(epsilon*np.cos(theta - thm) - 1.0)
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    m_in = (theta <= thm) & (y >= 0.0)
    m_out = (theta > thm) & (y >= 0.0)
    ax.plot([r[im]*np.cos(theta[im])], [r[im]*np.sin(theta[im])],
            's', ms=0.6, color='k')
    ax.plot(x[m_in], y[m_in], '-', color='gray', alpha=0.8, lw=0.5)
    ax.plot(x[m_out], y[m_out], '-', color='r', alpha=0.8, lw=0.5)
thm_grid = np.linspace(0.0, np.pi, 200)
rm = 2.0/(1.0 + np.cos(thm_grid))
xlocus = rm*np.cos(thm_grid)
ylocus = rm*np.sin(thm_grid)
ax.plot(xlocus, ylocus, '-', color='k', alpha=0.2, lw=3)

ax.plot([0.0], [0.0], '*', color='r')
ax.set(xlim=[-3.1, 3.9], ylim=[-0.1, 6.9],
       xlabel="$x / R_0$",
       ylabel="$y / R_0$")
sns.despine(trim=True)
fig.tight_layout()
fig.savefig(figfile)
print(figfile, end='')
