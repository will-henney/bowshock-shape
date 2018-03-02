INPUT=[["D theta", "Outer", "Ridge", "Inner"], [45, 1.775, 8.342, 3.312], [50, 1.775, 8.342, 2.771], [55, 2.153, 4.343, 3.163], [60, 2.186, 3.262, 3.163], [65, 2.395, 4.217, 2.692], [70, 2.478, 3.664, 2.911], [75, 2.749, 3.318, 2.721], [80, 2.749, 3.642, 2.858], ["Median", 2.291, 3.941, 2.885], ["MAD", 0.323, 0.513, 0.179]]
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
ax.legend(title="M42 000-400")
ax.axvspan(60.0, 80.0, color='0.9')
ax.set(ylim=[0.0, 9.0],
       xlabel=r"$\Delta \theta$, degrees",
       xticks=tab['D theta'],
       ylabel=r"Fitted planitude $\Pi'$")
sns.despine(trim=True)
fig.tight_layout()
fig.savefig(figfile)
print(figfile, end='')
