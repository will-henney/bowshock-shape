import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import equation6
import conic_parameters
sys.path.append('../conic-projection')
from conproj_utils import Conic
import crw_misc_utils


def polar_plot(r, th, ax, **kwargs):
    ax.plot(r*np.cos(th), r*np.sin(th), **kwargs)


try:
    beta = float(sys.argv[1])
except IndexError:
    sys.exit('Usage: {} BETA'.format(sys.argv[0]))


ntheta = 200

xmax = 1.5

R0 = np.sqrt(beta)/(1 + np.sqrt(beta))
shell = equation6.Shell(beta=beta)
theta = np.linspace(0.0, shell.th_infty, ntheta)
R, th1 = shell.radius(theta, method='brent', full=True)
R_approx = crw_misc_utils.radius(theta, crw_misc_utils.th1_approx(theta, beta))
m = R_approx > 0.0

th_tail = conic_parameters.theta_tail(beta, xi=None,
                                      f=conic_parameters.finf_CRW)
print('th1_infty =', np.degrees(shell.th1_infty), np.degrees(th_tail))

b_a = np.tan(th_tail)
x_tail = np.linspace(-xmax, xmax, 3)
y_tail = -b_a*(x_tail - 1.0)


shell2 = equation6.Shell(beta=beta, xi=1.0, innertype='anisotropic')
theta2 = np.linspace(0.0, shell2.th_infty, ntheta)
R2, th12 = shell2.radius(theta2, method='brent', full=True)
th_tail2 = conic_parameters.theta_tail(beta, xi=1.0)
print('th1_infty_2 =', np.degrees(shell2.th1_infty), np.degrees(th_tail2))

b_a2 = np.tan(th_tail2)
y_tail2 = -b_a2*(x_tail - 1.0)

figfilename = sys.argv[0].replace('.py', '-{:05d}.pdf').format(int(1e4*beta))
fig, (ax, axx) = plt.subplots(2, 1)
polar_plot(R, theta, ax)
polar_plot(R_approx[m], theta[m], ax, ls='None', marker='.', alpha=0.2)
polar_plot(R2, theta2, ax, lw=0.6)
ax.plot(x_tail, y_tail, '--')
ax.plot(x_tail, y_tail2, '--')

ax.set_xlim(-xmax, xmax)
ax.set_ylim(-0.2*xmax, 1.2*xmax)
ax.set_aspect('equal', adjustable='box')

axx.plot(shell.th_infty - theta, shell.th1_infty - th1)
axx.plot(shell2.th_infty - theta2, shell2.th1_infty - th12)


fig.set_size_inches(4.0, 6.0)
fig.savefig(figfilename)
print(figfilename)
