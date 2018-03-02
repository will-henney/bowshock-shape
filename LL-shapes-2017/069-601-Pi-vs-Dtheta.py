INPUT=[["D theta", "Outer", "Ridge", "Inner"], [45, 2.358, 4.354, 2.869], [50, 2.19, 4.354, 2.869], [55, 2.19, 3.045, 2.869], [60, 2.035, 2.619, 4.209], [65, 2.144, 2.496, 3.728], [70, 2.214, 2.384, 3.728], [75, 2.292, 2.355, 3.728], [80, 2.403, 2.275, 3.659], ["Median", 2.202, 2.558, 3.694], ["MAD", 0.074, 0.243, 0.275]]
import sys
import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
import seaborn as sns

tab = Table(rows=INPUT[1:9], names=INPUT[0])

figfile = sys.argv[0].replace('.py', '.pdf')

sns.set_color_codes()

fig, ax = plt.subplots(figsize=(4, 4))
arcs = "Inner", "Outer", "Ridge"
colors = 'mcr'
m = tab['D theta'] >= 60.0
for arc, color in zip(arcs, colors):
    mean = np.mean(tab[arc][m])
    sigma = np.std(tab[arc][m])
    label = fr"{arc}: ${mean:.1f} \pm {sigma:.1f}$"
    ax.plot(tab['D theta'], tab[arc], 'o', color=color, label=label)
ax.legend(title="M42 069-601")
ax.axvspan(60.0, 80.0, color='0.9')
ax.set(ylim=[0.0, 9.0],
       xlabel=r"$\Delta \theta$, degrees",
       xticks=tab['D theta'],
       ylabel=r"Fitted planitude, $\Pi'$")
sns.despine(trim=True)
fig.tight_layout()
fig.savefig(figfile)
print(figfile, end='')
