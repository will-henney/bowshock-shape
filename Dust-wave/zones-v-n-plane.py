import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from vec_root import chandrupatla

figname = sys.argv[0].replace('.py', '.pdf')

fig, axes = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(8, 4))
stardata = [
    [10.0, 0.63, 0.0026, axes[0]],
    [20.0, 5.5, 0.0446, axes[1]],
    [40.0, 22.0, 0.1682, axes[2]],
]

# Velocities in units of km/s (10 km/s -> 100 km/s)
vgrid = np.linspace(10.0, 100.0, 200)
# Densities in units of 1 pcc  (0.01 -> 1e5)
logngrid = np.linspace(-3.0, 6.0, 200)
# 2d versions of velocity and density grids
vv, nn = np.meshgrid(vgrid, 10**logngrid)

def rstar(v10, n, L4):
    """Characteristic radius in pc"""
    return 2.21*np.sqrt(L4/n)/v10

def taustar(v10, n, L4, kappa600=1.0):
    """Characteristic optical depth"""
    return 0.0089*kappa600*np.sqrt(L4*n)/v10

def xfunc(x, ts, eta):
    """Function to be zeroed to find x"""
    return x**2 - (1.0 - np.exp(-2*ts*x)) - eta

R0s = (0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0)
for M, L4, eta, ax in stardata:
    Rs = rstar(vv/10, nn, L4)
    ts = taustar(vv/10, nn, L4)
    a, b = 0.0, 2*np.sqrt(1.0 + eta)
    x = chandrupatla(xfunc, a, b, args=(ts, eta))
    R0 = x*Rs
    tau = 2*x*ts
    ax.contour(vv, nn, tau, (eta/3, eta, 3*eta), colors='r')
    ax.contour(vv, nn, tau, (1.0, 3.0), colors='m')
    ax.contour(vv, nn, R0, R0s, colors='b')
    ax.set(yscale='log')

axes[0].set(xlabel=r"$v$, km/s", ylabel=r"$n$, pcc")
fig.tight_layout()
fig.savefig(figname)

print(figname, end='')
