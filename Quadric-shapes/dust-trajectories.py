import sys
import numpy as np
from matplotlib import pyplot as plt

figfile = sys.argv[0].replace('.py', '.pdf')
NTH = 501
fig, ax = plt.subplots(figsize=(5, 5))
for thm in np.linspace(0.0, 1.5*np.pi/2, 25):
    rm = 2.0/(1.0 + np.cos(thm))
    epsilon = 1./np.cos(thm)
    theta = np.linspace(0.01, 2*thm - 0.01, NTH)
    r = 0.5*(epsilon**2 - 1)/(epsilon*np.cos(theta - thm) - 1.0)
    ax.plot(r*np.cos(theta), r*np.sin(theta), '-', color='k', lw=0.5)
    ax.plot([rm*np.cos(thm)], [rm*np.sin(thm)], '.', color='k')

ax.plot([0.0], [0.0], '*', color='r')
ax.set(xlim=[-1.5, 2.5], ylim=[-0.2, 3.8])
fig.savefig(figfile)
print(figfile, end='')
