import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from dust_bfield_ode import streamline

figfile = sys.argv[0].replace('.py', '.pdf')

# Impact parameter
Y0 = 0.35
Z0 = 0.1
thB_degrees = -28.2
stream = streamline(Y0=Y0, Z0=Z0, thB=np.radians(thB_degrees), tstop=150, X0=10., n=2001)
sns.set_style('white')
sns.set_color_codes()
fig, (ax, axp) = plt.subplots(2, 1, figsize=(4, 6))
ax.plot(stream['t'], stream['u'], label='$U$')
ax.plot(stream['t'], stream['v'], label='$V$')
ax.plot(stream['t'], stream['x'], label='$X$')
ax.plot(stream['t'], stream['y'], label='$Y$')
ax.axhspan(0.0, 1.0, color='k', alpha=0.1)
label = fr"$\theta_B = {thB_degrees:.1f}^\circ$, "
label += f"$y_0 = {Y0:.2f}$, $z_0 = {Z0:.2f}$"
ax.legend(title=label, ncol=2)
ax.set(xlabel='Time', ylim=[-2.5, 4.5], xlim=[0, 70])

axp.plot(stream['x'], stream['u'], label='$(X, U)$')
axp.plot(stream['y'], stream['v'], label='$(Y, V)$')
axp.axhline(0, color='k', lw=0.5)
axp.legend(title='Phase space')
axp.set(xlabel='$X$, $Y$', ylabel='$U$, $V$',
        xlim=[-7, 9], ylim=[-1.1, 1.1])

sns.despine(trim=True)
fig.tight_layout()
fig.savefig(figfile)
print(figfile, end='')
