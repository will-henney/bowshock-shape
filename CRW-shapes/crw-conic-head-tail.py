import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import equation6
import conic_parameters
sys.path.append('../conic-projection')
from conproj_utils import Conic

def R_from_theta_CRW(theta, beta):
    shell = equation6.Shell(innertype='isotropic', beta=beta)
    R = shell.radius(theta)
    return R/shell.R0

betagrid = [1e-4, 0.001, 0.01, 0.02, 0.03, 0.1, 0.5]
nbeta = len(betagrid)

ntheta = 100
theta = np.linspace(0.0, np.pi, ntheta)

figfilename = sys.argv[0].replace('.py', '.pdf')
fig, axes = plt.subplots(1, nbeta, sharex=True, sharey=True)

ytop = 6.5
baseline = 0.7
for i, beta in enumerate(betagrid):
    ax = axes[i]

    # Vanilla CRW solution
    R_crw = R_from_theta_CRW(theta, beta)
    x_crw = R_crw*np.cos(theta)
    y_crw = R_crw*np.sin(theta)

    ht = conic_parameters.HeadTail(beta)

    x_head = ht.x_head(ht.t_h)
    y_head = ht.y_head(ht.t_h)

    x_tail = ht.x_tail(ht.t_t)
    y_tail = ht.y_tail(ht.t_t)

    ax.plot(x_crw, y_crw)
    ax.plot(x_head, y_head, '--')
    ax.plot(x_tail, y_tail, '--')
    ax.plot(0.0, 0.0, 'o', color='k')
    ax.axhline(ls=':')
    ax.axvline(ls=':')
    ax.set_aspect('equal', adjustable='box-forced')
    ax.text(2.0, ytop,
            r'$\beta = {:.4f}$ '.format(beta), ha='right')
    ax.text(2.0, ytop - baseline,
            r'$x_m = {:.2f}$ '.format(ht.x_m), ha='right')
    ax.text(2.0, ytop - 2*baseline,
            r'$x_t = {:.2f}$ '.format(ht.x0_t), ha='right')
    ax.text(2.0, ytop - 3*baseline,
            r'$x_h = {:.2f}$ '.format(ht.x0_h), ha='right')
    ax.text(2.0, ytop - 4*baseline,
            r"$\phi' / \phi = {:.2f}$ ".format(ht.phi1_over_phi), ha='right')
    ax.text(2.0, ytop - 5*baseline,
            r"$\theta_h = {:.2f}$ ".format(np.degrees(ht.theta_h)), ha='right')
    ax.text(2.0, ytop - 6*baseline,
            r"$\theta_t = {:.2f}$ ".format(np.degrees(ht.theta_t)), ha='right')


axes[0].set_xlim(-5.0, 2.1)
axes[0].set_ylim(-0.1, 7.0)
axes[0].set_xlabel(r'$x / r_{0}$')
axes[0].set_ylabel(r'$y / r_{0}$')

fig.set_size_inches(3*nbeta, 3)
fig.tight_layout()
fig.savefig(figfilename)
print(figfilename)
