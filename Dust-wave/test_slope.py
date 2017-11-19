import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from bow_projection import (omega, paraboloid_R_theta, alpha,
                            wilkinoid_R_theta, cantoid_R_theta,
                            Spline_R_theta_from_function)

figfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('ticks')
fig, ax = plt.subplots()

th = np.linspace(0, np.pi, 10001)
th_dg = np.degrees(th)

for label, func, pars, ngrid, s in [
        ["paraboloid", paraboloid_R_theta, (), 1001, 0],
        ["Wilkinoid", wilkinoid_R_theta, (), 1001, 0],
        [r"Cantoid $\beta = 0.001$", cantoid_R_theta, (0.001,), 1001, 0],
        [r"Cantoid $\beta = 0.01$", cantoid_R_theta, (0.01,), 1001, 0],
        [r"Cantoid $\beta = 0.1$", cantoid_R_theta, (0.1,), 1001, 0],
]:
    spline_func = Spline_R_theta_from_function(
        ngrid=ngrid, smooth=s, shape_func=func, shape_func_pars=pars)
    alpha_slope = np.degrees(alpha(th, func, *pars))
    alpha_slope_spline = np.degrees(alpha(th, spline_func))
    ax.plot(th_dg, (90 - alpha_slope)/th_dg, 
            color='b', alpha=0.2, lw=2, label='_nolabel_')
    ax.plot(th_dg, (90 - alpha_slope_spline)/th_dg,
            lw=0.8, label=label)

ax.legend(title=r"Spline approximations")
ax.set(
    xlabel=r"Polar angle: $\theta$, degrees",
    ylabel=r"Slope angle, $(90 - \alpha) / \theta$",
    xlim=[0, 180],
    ylim=[0.0, 1.1],
    xticks=[0, 30, 60, 90, 120, 150, 180],
)
sns.despine()
fig.tight_layout()
fig.savefig(figfile)
print(figfile, end='')
