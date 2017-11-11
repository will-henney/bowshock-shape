import sys
import numpy as np
from scipy.optimize import fsolve
from matplotlib import pyplot as plt
import matplotlib.ticker
from matplotlib.patches import Ellipse
import seaborn as sns
import conic_parameters
from equation6 import Shell

# Vectorized version of this function since the one in
# conic_parameters is scalar-only
def theta_tail(beta, xi, f=conic_parameters.finf, th_init=np.radians(91.0)):
    thinf = np.empty_like(beta)
    for i, b in enumerate(beta):
        thinf[i] = fsolve(f, th_init, args=(b,xi))
    return np.pi - thinf


plotfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('ticks')
sns.set_color_codes(palette='deep')
fig, ax = plt.subplots(figsize=(4, 5))

ks = [0.0, 0.5, 3.0, 8.0]
colors = 'krmb'
beta = np.logspace(-4.0,-0.01,1000)

xmin, xmax = 0.0001, 1.0
ymin, ymax = -90.0, 90.0
labelsize = "small"
legend_size = "medium"

#ax.fill_between([xmin, xmax], [ymin, ymin], [0, 0], color='y', alpha=0.05)
ax.fill_between([xmin, xmax], [0, 0], [45, 45], color='k', alpha=0.2)
ax.fill_between([xmin, xmax], [45, 45], [ymax, ymax], color='k', alpha=0.1)
wbox=dict(facecolor='white', alpha=0.9, ec='none', pad=1.0)
#ax.text(0.04, 47.0, 'Oblate spheroid', bbox=wbox, fontsize=labelsize)
ax.text(0.1, 70.0, 'Oblate\nspheroids',
        ha='center', bbox=wbox, fontsize=labelsize)
ax.text(0.1, 12.0, 'Prolate\nspheroids',
        ha='center', bbox=wbox, fontsize=labelsize)
ax.text(0.1, -12.0, 'Hyperboloids',
        ha='center', bbox=wbox, fontsize=labelsize)
ax.axhline(0.0, color='k', lw=0.5)
ax.text(0.1, 45.0, 'Sphere', va='center',
        ha='center', bbox=wbox, fontsize=labelsize)
ax.text(0.1, 0.0, 'Paraboloid', va='center',
        ha='center', bbox=wbox, fontsize=labelsize)
for k, color in zip(ks,colors):
    xi = 2.0 / (2.0 + k)
    thc = np.degrees(conic_parameters.theta_c(beta, xi))
    thct = np.degrees(-theta_tail(beta, xi))
    ax.plot(beta, thc, linestyle="-", color=color, label=fr"$k={k}$")
    ax.plot(beta, thct, linestyle="--", color=color, label=None)

thct_crw = np.degrees(-theta_tail(beta, None, f=conic_parameters.finf_CRW))
ax.plot(beta, thct_crw, linestyle="-.", color="k", label="_nolabel_")

i0 = len(beta)//2
ax.annotate("Cantoids", xy=(beta[i0], thct_crw[i0]),
            xytext=(-30, -30), textcoords='offset points',
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))

xy_head = (4e-4, 33.0)
xy_tail = (4e-4, -5.0)
ellipse_props = dict(ec='k', fc='w', zorder=100, alpha=0.5)
ellipse_head = Ellipse(xy_head, 1e-4, 20, **ellipse_props)
ellipse_tail = Ellipse(xy_tail, 1e-4, 20, **ellipse_props)
ax.add_patch(ellipse_head)
ax.add_patch(ellipse_tail)
ax.annotate(r"$\theta_{\mathcal{Q}}$ from $\Pi$, $\Lambda$", xy=xy_head,
            xytext=(15, 30), textcoords='offset points',
            arrowprops=dict(arrowstyle="simple,tail_width=0.1",
                            fc="0.2", ec="none",
                            patchB=ellipse_head, shrinkB=17,
                            connectionstyle="angle3,angleA=-10,angleB=-90"),
)
ax.annotate(r"$\theta_{\mathcal{Q}}$ from $\theta_\infty$", xy=xy_tail,
            xytext=(10, -40), textcoords='offset points',
            arrowprops=dict(arrowstyle="simple,tail_width=0.1",
                            fc="0.2", ec="none",
                            patchB=ellipse_tail, shrinkB=17,
                            connectionstyle="angle3,angleA=0,angleB=-100"),
)

leg = ax.legend(loc='lower left', fontsize=legend_size, ncol=2,
                title='Anisotropy index')
leg.get_title().set_fontsize(legend_size)
ax.set(
    xlabel=r"Axial momentum ratio: $\beta$",
    ylabel=r"Equivalent quadric angle, $\theta_{\mathcal{Q}}$",
    xscale='log',
    xlim=[xmin, xmax],
    ylim=[ymin, ymax],
    yticks=[-90, -60, -30, 0, 30, 60, 90],
)

sns.despine()
fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
