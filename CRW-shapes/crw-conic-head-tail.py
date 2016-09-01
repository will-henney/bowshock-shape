import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import equation6
import conic_parameters
sys.path.append('../conic-projection')
from conproj_utils import Conic

try:
    XI_PARAMETER = float(sys.argv[1])
    suffix = '-xi{:03d}'.format(int(100*XI_PARAMETER))
except IndexError:
    XI_PARAMETER = None
    suffix = ''


def R_from_theta_CRW(theta, beta):
    shell = equation6.Shell(innertype='isotropic', beta=beta)
    R = shell.radius(theta)
    return R/shell.R0

betagrid = [3e-5, 1e-4, 0.001, 0.003, 0.01, 0.03, 0.1, 0.5]
nbeta = len(betagrid)

ntheta = 100
theta = np.linspace(0.0, np.pi, ntheta)

figfilename = sys.argv[0].replace('.py', suffix + '.pdf')

sns.set_style('whitegrid')
sns.set_color_codes('dark')

NROWS = 2
fig, axes = plt.subplots(NROWS, nbeta//NROWS, sharex=True, sharey=True)

ytop = 10.5
xright = 4.0
baseline = 0.7
whitebox = {'edgecolor': 'none', 'facecolor': 'white',
            'alpha': 0.7, 'boxstyle': 'round,pad=0.1'}

xmin, xmax = -7.0, 4.1
ymin, ymax = -0.1, 11.0

# x-data for tail asymptote
xa = np.linspace(xmin, xmax, 2)

for i, beta in enumerate(betagrid):
    ax = axes.flat[i]

    # Vanilla CRW solution
    R_crw = R_from_theta_CRW(theta, beta)
    x_crw = R_crw*np.cos(theta)
    y_crw = R_crw*np.sin(theta)

    ht = conic_parameters.HeadTail(beta, xi=XI_PARAMETER, xmin=-1.0)

    x_head = ht.x_head(ht.t_h)
    y_head = ht.y_head(ht.t_h)

    x_tail = ht.x_tail(ht.t_t)
    y_tail = ht.y_tail(ht.t_t)

    # asymptote to tail
    ya = (ht.x0_t - xa)*np.tan(ht.theta_t)
    ax.plot(xa, ya, lw=0.3, color='orange')

    ax.plot(x_crw, y_crw, lw=3, color='y', alpha=0.5)
    ax.plot(x_head, y_head, '--', color='g')
    ax.plot(x_tail, y_tail, '-', lw=1.0, color='r')
    ax.plot(0.0, 0.0, 'o', color='k')
    ax.plot(ht.x0_t - ht.a_t, 0.0, '.', color='k')
    ax.axhline(ls=':')
    ax.axvline(ls=':')
    ax.axvline(x=ht.x_m, lw=1.0, ls='--', color='black', alpha=0.3)
    ax.set_aspect('equal', adjustable='box-forced')
    text = r'$\beta = {:.5f}$ '.format(beta)
    text += '\n' + r'$x_m = {:.2f}$ '.format(ht.x_m)
    text += '\n' + r'$x_t = {:.2f}$ '.format(ht.x0_t)
    text += '\n' + r'$x_t - a_t = {:.2f}$ '.format(ht.x0_t - ht.a_t)
    text += '\n' + r'$x_h = {:.2f}$ '.format(ht.x0_h)
    text += '\n' + r"$J = {:.2f}$ ".format(ht.J)
    # text += '\n' + r"$K = {:.2f}$ ".format(ht.K)
    text += '\n' + r"$\theta_h = {:.2f}$ ".format(np.degrees(ht.theta_h))
    text += '\n' + r"$\theta_t = {:.2f}$ ".format(np.degrees(ht.theta_t))
    ax.text(xright, ytop, text, ha='right', va='top', bbox=whitebox, fontsize='small')


axes[1, 0].set_xlim(xmin, xmax)
axes[1, 0].set_ylim(ymin, ymax)
axes[1, 0].set_xlabel(r'$x / r_{0}$')
axes[1, 0].set_ylabel(r'$y / r_{0}$')

fig.set_size_inches(3*nbeta/NROWS, 3*NROWS)
fig.tight_layout()
fig.savefig(figfilename)
print(figfilename)
