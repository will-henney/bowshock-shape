import sys
import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt
import seaborn as sns

def dydt(y, t, alpha):
    """Derivatives for ODE: x'' = 0.5 (x^{-2} - alpha^2 (x' + 1))"""
    x, u = y
    dxdt = u
    dudt = 0.5*(x**(-2) - alpha**2 * (u + 1.0))
    return [dxdt, dudt]

figfile = sys.argv[0].replace('.py', '.pdf')

# Initial conditions
y0 = [10.0, -1.0]

# Coupling parameter
alpha = 1.0/2.0

# Time grid
t = np.linspace(0.0, 60.0, 501)
soln = odeint(dydt, y0, t, args=(alpha,))

# Slippage velocity
w = 1.0 + soln[:, 1]
# Drift velocity
wdrift = 1.0 / alpha / soln[:, 0]

sns.set_style('ticks')
sns.set_color_codes('dark')
fig, (ax, axp) = plt.subplots(2, 1, figsize=(4, 6))
ax.plot(t, soln[:, 0], label='$x$')
ax.plot(t, w, label='$w$')
ax.plot(t, wdrift, ls='--', label='$w_\mathrm{drift}$')

ax.axhline(1.0/alpha, ls=':', color='k', lw=0.8)
ax.axhspan(0.0, 1.0, color='k', alpha=0.1)
ax.legend(title=r"$\alpha_\mathrm{drag} = 0.5$")
ax.set(xlabel='Time', ylim=[None, 4])
t2 = np.linspace(0.0, 20.0, 201)
soln2 = odeint(dydt, y0, t2, args=(2.0,))
soln0 = odeint(dydt, y0, t2, args=(0.0,))

axp.plot(soln0[:, 0], soln0[:, 1], label=r"$\alpha_\mathrm{drag} = 0$")
axp.plot(soln[:, 0], soln[:, 1], label=r"$\alpha_\mathrm{drag} = 0.5$")
axp.plot(soln2[:, 0], soln2[:, 1], label=r"$\alpha_\mathrm{drag} = 2$")
axp.axhline(0, xmax=0.55, color='k', lw=0.5)
axp.legend(title='Phase space')
axp.set(xlabel='$x$', ylabel='$u$',
        xlim=[0, 7], ylim=[-1, 1], yticks=[-1.0, -0.5, 0., 0.5, 1.0])

sns.despine()
fig.tight_layout()
fig.text(0.02, 0.97, '(a)')
fig.text(0.02, 0.5, '(b)')
fig.savefig(figfile)
print(figfile, end='')
