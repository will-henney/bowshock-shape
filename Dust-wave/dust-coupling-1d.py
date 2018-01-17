import sys
import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt
import seaborn as sns

def dydt(y, t, alpha, sound=0.5):
    """Derivatives for ODE: x'' = 0.5 (x^{-2} - alpha^2 (x' + 1))"""
    x, u = y
    dxdt = u
    dudt = 0.5*(x**(-2)
                - alpha**2 * (u + 1.0)*np.sqrt((1.65*sound)**2
                                               + (u + 1.0)**2))
    return [dxdt, dudt]

figfile = sys.argv[0].replace('.py', '.pdf')

# Initial conditions
y0 = [10.0, -1.0]

# Coupling parameter
alpha = 0.4

# Time grid
t = np.linspace(0.0, 60.0, 501)
soln = odeint(dydt, y0, t, args=(alpha,))
t0 = t[np.argmin(soln[:, 0])]

# Slippage velocity
w = 1.0 + soln[:, 1]
# Drift velocity
wdrift = 1.0 / alpha / soln[:, 0]

sns.set_style('ticks')
sns.set_color_codes('dark')
fig, (ax, axp) = plt.subplots(2, 1, figsize=(4, 6))
ax.plot(t - t0, soln[:, 0], label='$R/R_{0}$')
ax.plot(t - t0, w, label='$w / v_{\infty}$')
ax.plot(t - t0, wdrift, ls='--', label='$w_\mathrm{drift} / v_{\infty}$')

ax.axhline(1.0/alpha, ls=':', color='k', lw=0.8)
ax.axhspan(0.0, 1.0, color='k', alpha=0.1)
ax.legend(title=r"$\alpha_\mathrm{drag} = 0.5$")
ax.set(xlabel=r'Time / $(R_{0} / v_{\infty})$', ylim=[-0.3, 4.3])
t2 = np.linspace(0.0, 20.0, 201)
soln1 = odeint(dydt, y0, t2, args=(1.0,))
soln2 = odeint(dydt, y0, t2, args=(2.0,))
soln0 = odeint(dydt, y0, t2, args=(0.0,))

axp.plot(soln0[:, 0], soln0[:, 1], label=r"$\alpha_\mathrm{drag} = 0$")
axp.plot(soln[:, 0], soln[:, 1], label=r"$\alpha_\mathrm{drag} = 0.5$")
axp.plot(soln1[:, 0], soln1[:, 1], label=r"$\alpha_\mathrm{drag} = 1$")
axp.plot(soln2[:, 0], soln2[:, 1], label=r"$\alpha_\mathrm{drag} = 2$")
axp.axhline(0, xmax=0.55, color='k', lw=0.5)
axp.legend(title='Phase space\n  trajectories')
axp.set(xlabel='$R/R_{0}$', ylabel='$v / v_{\infty}$',
        xlim=[-0.35, 6.9], ylim=[-1.1, 1.1],
        xticks=range(7),
        yticks=[-1.0, -0.5, 0., 0.5, 1.0])

sns.despine(trim=True)
fig.tight_layout()
fig.text(0.02, 0.97, '(a)')
fig.text(0.02, 0.5, '(b)')
fig.savefig(figfile)
print(figfile, end='')
