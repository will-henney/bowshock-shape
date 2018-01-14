import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from dust_blorentz_ode import streamline

figfile = sys.argv[0].replace('.py', '.pdf')

# Impact parameter
Y0 = 0.00
Z0 = 0.00
thB_degrees = 25
MACH_ALFVEN = 4.0
stream = streamline(Y0=Y0, Z0=Z0, thB=np.radians(thB_degrees), tstop=150, X0=10., n=2001, LFAC=10, ALPHA_DRAG=0.5, V_TURB_0=1.0/MACH_ALFVEN)
sns.set_style('white')
sns.set_color_codes()
fig, (ax, axp) = plt.subplots(2, 1, figsize=(4, 6))
ax.plot(stream['t'], stream['u'], label='$U$')
ax.plot(stream['t'], stream['v'], label='$V$')
ax.plot(stream['t'], stream['w'], label='$W$')
ax.plot(stream['t'], stream['x'], label='$X$')
ax.plot(stream['t'], stream['y'], label='$Y$')
ax.plot(stream['t'], stream['z'], label='$Z$')
ax.axhspan(0.0, 1.0, color='k', alpha=0.1)
label = fr"$\theta_B = {thB_degrees:.1f}^\circ$, "
label += f"$y_0 = {Y0:.2f}$, $z_0 = {Z0:.2f}$"
ax.legend(title=label, ncol=2)
ax.set(xlabel='Time', ylim=[-2, 3], xlim=[0, 50])

axp.plot(stream['x'], stream['u'], label='$(X, U)$')
axp.plot(stream['y'], stream['v'], label='$(Y, V)$')
axp.plot(stream['z'], stream['w'], label='$(Z, W)$')
axp.axhline(0, color='k', lw=0.5)
axp.legend(title='Phase space')
axp.set(xlabel='$X$, $Y$, $Z$', ylabel='$U$, $V$, $W$',
        xlim=[-3, 3], ylim=[-1.5, 1.5])

sns.despine(trim=True)
fig.tight_layout()
fig.savefig(figfile)
print(figfile, end='')
