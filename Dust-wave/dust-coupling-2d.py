import sys
from matplotlib import pyplot as plt
import seaborn as sns
from dust_couple_ode import streamline

figfile = sys.argv[0].replace('.py', '.pdf')

# Impact parameter
b = 0.001
stream = streamline(Y0=b, alpha=1.0/2.0)
sns.set_style('white')
sns.set_color_codes()
fig, (ax, axp) = plt.subplots(2, 1, figsize=(4, 6))
ax.plot(stream['t'], stream['u'], label='$u$')
ax.plot(stream['t'], stream['v'], label='$v$')
ax.plot(stream['t'], stream['x'], label='$x$')
ax.plot(stream['t'], stream['y'], label='$y$')
ax.axhspan(0.0, 1.0, color='k', alpha=0.1)
ax.legend(title=(fr"$\alpha = {stream['alpha']:.3f}$, "
                 f"$b = {b:.3f}$"), ncol=2)
ax.set(xlabel='Time', ylim=[-3, 5])

axp.plot(stream['x'], stream['u'], label='$(x, u)$')
axp.plot(stream['y'], stream['v'], label='$(y, v)$')
axp.axhline(0, color='k', lw=0.5)
axp.legend(title='Phase space')
axp.set(xlabel='$x$, $y$', ylabel='$u$, $v$',
        xlim=[-3, 5], ylim=[-1.1, 1.1])

sns.despine()
fig.tight_layout()
fig.savefig(figfile)
print(figfile, end='')
