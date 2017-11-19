import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import bow_projection as bp
import standing_wave

figfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('ticks')
fig, ax = plt.subplots()

th = np.linspace(-np.pi, np.pi, 1001)
th_dg = np.degrees(th)

amplitude, wavenumber = 0.001, 20.0
shape = standing_wave.StandingWave(
    bp.wilkinoid_R_theta, amplitude=amplitude, wavenumber=wavenumber)
phases = [0.0, 0.25, 0.5]
for phase in phases:
    shape.phase = phase
    alpha_slope = np.degrees(bp.alpha(th, shape))
    ax.plot(th_dg, (90.0 - alpha_slope)/th_dg,
            lw=0.8, label=f"phase {phase:.2f}")

ax.legend(title=rf"$A = {amplitude:.1f}$, $N = {wavenumber:.1f}$")
ax.set(
    xlabel=r"Polar angle: $\theta$, degrees",
    ylabel=r"Slope angle, $(90 - \alpha)/\theta$",
    xlim=[0, 180],
    ylim=[0, 1.1],
    xticks=[0, 30, 60, 90, 120, 150, 180],
)
sns.despine()
fig.tight_layout()
fig.savefig(figfile)
print(figfile, end='')
