import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from bow_projection import sin_phi_t, paraboloid_R_theta, theta_infinity, theta_0_90

figfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('ticks')
fig, ax = plt.subplots(figsize=(5, 5))

inclinations = [0, 15, 30, 45, 60, 75]
colors = sns.color_palette(n_colors=len(inclinations))
th_inf = theta_infinity(paraboloid_R_theta)
for inc_dg, color in zip(inclinations, colors):
    inc = np.radians(inc_dg)
    th0, th90 = theta_0_90(inc, paraboloid_R_theta)
    th = np.linspace(th0, th_inf, 501)
    th_dg = np.degrees(th)
    print("theta:", th_dg, file=sys.stderr)
    sphit = sin_phi_t(th, inc, paraboloid_R_theta)
    print("sphit:", sphit, file=sys.stderr)
    phit_dg = np.degrees(np.arcsin(sphit))
    print("phi:", phit_dg, file=sys.stderr)
    ax.plot(th_dg, phit_dg, label=f"inc = {inc_dg:d}", color=color)

ax.legend(title="paraboloid")
ax.set(
    xlabel=r"$\theta$",
    ylabel=r"$\phi_t$",
    xlim=[0, 180],
    ylim=[-90, 90],
    xticks=[0, 30, 60, 90, 120, 150, 180],
    yticks=[-90, -60, -30, 0, 30, 60, 90],
)
sns.despine()
fig.tight_layout()
fig.savefig(figfile)
print(figfile, end='')
