import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from bow_projection import sin_phi_t, paraboloid_R_theta

figfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('ticks')
fig, ax = plt.subplots(figsize=(5, 5))

th = np.linspace(0.0, np.pi, 501)
th_dg = np.degrees(th)
inclinations = [0, 15, 30, 45, 60, 75]
colors = sns.color_palette(n_colors=len(inclinations))
for inc_dg, color in zip(inclinations, colors):
    inc = np.radians(inc_dg)
    sphit = sin_phi_t(th, inc, paraboloid_R_theta)
    phit_dg = np.degrees(np.arcsin(sphit))
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
