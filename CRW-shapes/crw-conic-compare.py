import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import equation6
import conic_parameters
sys.path.append('../conic-projection')
from conproj_utils import Conic


def R_from_theta(theta, beta, xi):
    shell = equation6.Shell(innertype='anisotropic', beta=beta, xi=xi)
    R = shell.radius(theta)
    R0 = np.sqrt(beta)/(1 + np.sqrt(beta))
    return R/R0


xigrid = [1.0, 0.8, 0.4, 0.1]
nxi = len(xigrid)
betagrid = [0.001, 0.01, 0.1]
nbeta = len(betagrid)

ntheta = 100
theta = np.linspace(0.0, np.pi, ntheta)

figfilename = sys.argv[0].replace('.py', '.pdf')
fig, axes = plt.subplots(nxi, nbeta, sharex=True, sharey=True)

for j, xi in enumerate(xigrid):
    for i, beta in enumerate(betagrid):
        ax = axes[j, i]

        # Geberalized CRW solution
        R_crw = R_from_theta(theta, beta, xi)
        x_crw = R_crw*np.cos(theta)
        y_crw = R_crw*np.sin(theta)

        # Matched conic parameters
        A = conic_parameters.A(beta, xi)
        th_conic = np.degrees(conic_parameters.theta_c(beta, xi))
        c = Conic(A=A, th_conic=th_conic)
        t = c.make_t_array()
        x_con = c.x(t)
        y_con = c.y(t)

        # Compare the two curves
        ax.plot(x_crw, y_crw)
        ax.plot(x_con, y_con, '--')
        ax.plot(0.0, 0.0, 'o', color='k')
        ax.axhline(ls=':')
        ax.axvline(ls=':')
        ax.text(1.5, 2.5, r'$\beta = {:.3f}$ '.format(beta), ha='right')
        ax.text(1.5, 2.0, r'$\xi = {:.1f}$ '.format(xi), ha='right')

axes[-1, 0].set_xlim(-1.5, 1.5)
axes[-1, 0].set_ylim(-0.1, 3.0)
axes[-1, 0].set_xlabel(r'$x / r_{0}$')
axes[-1, 0].set_ylabel(r'$y / r_{0}$')

fig.set_size_inches(9, 12)
fig.tight_layout()
fig.savefig(figfilename)
print(figfilename)
