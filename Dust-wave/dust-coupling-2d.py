import sys
from matplotlib import pyplot as plt
import seaborn as sns
from dust_couple_ode import streamline

figfile = sys.argv[0].replace('.py', '.pdf')

# Impact parameter
b = 0.001
stream = streamline(Y0=b, alpha=1.0/2.0, tstop=70, X0=20., n=501, mu=0.01)
sns.set_style('white')
sns.set_color_codes()
fig, (ax, axp) = plt.subplots(2, 1, figsize=(4, 6))
ax.plot(stream['t'], stream['u'], label='$U$')
ax.plot(stream['t'], stream['v'], label='$V$')
ax.plot(stream['t'], stream['x'], label='$X$')
ax.plot(stream['t'], stream['y'], label='$Y$')
ax.axhspan(0.0, 1.0, color='k', alpha=0.1)
label = fr"$\alpha = {stream['alpha']:.2f}$, "
if stream['mu'] is not None:
    label += fr"$\mu = {stream['mu']:.2f}$, "
label += f"$b = {b:.3f}$"
ax.legend(title=label, ncol=2)
ax.set(xlabel='Time', ylim=[-6, 8])

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
