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
t = np.linspace(0.0, 60.0, 201)

soln = odeint(dydt, y0, t, args=(alpha,))

sns.set_style('white')
sns.set_color_codes('dark')
fig, (ax, axp) = plt.subplots(2, 1, figsize=(4, 6))
ax.plot(t, soln[:, 1], label='$u$')
ax.plot(t, soln[:, 0], label='$x$')
ax.axhspan(0.0, 1.0, color='k', alpha=0.1)
ax.legend()
ax.set(xlabel='Time', ylim=[None, 4])

axp.plot(soln[:, 0], soln[:, 1])
axp.axhline(0, color='k', lw=0.5)
axp.set(xlabel='$x$', ylabel='$u$',
        xlim=[None, 5], ylim=[-1, 1])

sns.despine()
fig.tight_layout()
fig.savefig(figfile)
print(figfile, end='')
