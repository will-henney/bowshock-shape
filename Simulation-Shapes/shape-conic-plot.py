import sys
import json
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

try:
    prefix = sys.argv[1]
except:
    print(f"Usage: {sys.argv[0]} PREFIX")

jfile = f'{prefix}-arcdata.json'
data = json.load(open(jfile))
plotfile = sys.argv[0].replace('.py', f'-{prefix}.pdf')

R0 = np.array(data['outer']['R0'])
R = np.array(data['outer']['R'])
th = np.radians(data['outer']['theta'])

sns.set_style('ticks')
fig, ax = plt.subplots()
ax.plot(np.cos(th), R0/R, 'o')
ax.plot([-1, 1], [0, 1], '-', c='k', lw=0.5)
ax.set(
    xlim=[-1.05, 1.05],
    ylim=[-0.05, 1.05],
    xlabel=r"$\cos \,\theta$",
    ylabel="$R_{0} / R$",
)
sns.despine(trim=True)
fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
