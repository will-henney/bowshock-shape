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


def omega(R, theta):
    """(1/R) (d R/d theta) by finite differences"""
    dR_dtheta = np.pad(np.diff(R)/np.diff(theta), (0, 1), 'edge')
    return dR_dtheta/R


def alpha(R, theta):
    """Angle of tangent to curve R(theta)"""
    t = np.tan(theta)
    om = omega(R, theta)
    tana = (1 + om*t)/(t - om)
    return np.arctan(tana)


try:
    beta = float(sys.argv[1])
except IndexError:
    sys.exit('Usage: {} BETA'.format(sys.argv[0]))

try:
    fix_K = float(sys.argv[2])
except IndexError:
    fix_K = None


ntheta = 400
nearly = 1.0 # - 1e-5

xmax = 1.5

R0 = np.sqrt(beta)/(1 + np.sqrt(beta))
shell = equation6.Shell(beta=beta)
theta = np.linspace(0.0, nearly*shell.th_infty, ntheta)
R, th1 = shell.radius(theta, method='brent', full=True)
alph = alpha(R, theta)
R_approx = crw_misc_utils.radius(theta, crw_misc_utils.th1_approx(theta, beta))
m = R_approx > 0.0

th_tail = conic_parameters.theta_tail(beta, xi=None,
                                      f=conic_parameters.finf_CRW)
print('th1_infty =', np.degrees(shell.th1_infty), np.degrees(th1[-4:]))
print('alpha_infty =', np.degrees(alph[-4:]))

# Gradient: d phi_1 / d phi
grad = np.diff(shell.th1_infty - th1) / np.diff(shell.th_infty - theta)
# Theoretical estimate:
grad0 = beta*(np.pi / (shell.th1_infty
                       - np.sin(shell.th1_infty)*np.cos(shell.th1_infty)) - 1)
print('gradient:', grad0, grad[-4:])

b_a = np.tan(th_tail)
x_tail = np.linspace(-xmax, xmax, 3)
y_tail = -b_a*(x_tail - 1.0)


# 30 Aug 2016 - add in the attempted quadratic fit to phi_1 vs phi
ht = conic_parameters.HeadTail(beta)
print('Original tail parameters:')
print('beta = {:.4f}, tau = {:.2f}, J = {:.2f}, K = {:.2f}'.format(beta, ht.tau_t, ht.J, ht.K))
# ht.K *= ht.tau_t**2
if fix_K is not None: 
    ht.K = fix_K
    print('Corrected K = {:.2f}'.format(ht.K))

def fquad(phi, J=ht.J, K=ht.K):
    return J*phi + K*phi**2

shell2 = equation6.Shell(beta=beta, xi=1.0, innertype='anisotropic')
theta2 = np.linspace(0.0, nearly*shell2.th_infty, ntheta)
R2, th12 = shell2.radius(theta2, method='brent', full=True)
alph2 = alpha(R2, theta2)
th_tail2 = conic_parameters.theta_tail(beta, xi=1.0)
print('th1_infty_2 =', np.degrees(shell2.th1_infty), np.degrees(th12[-4:]))
print('alpha_infty_2 =', np.degrees(alph2[-4:]))

b_a2 = np.tan(th_tail2)
y_tail2 = -b_a2*(x_tail - 1.0)

figfilename = sys.argv[0].replace('.py', '-{:05d}.pdf').format(int(1e4*beta))
fig, (ax, axx, axxx) = plt.subplots(3, 1)
polar_plot(R, theta, ax)
polar_plot(R_approx[m], theta[m], ax, ls='None', marker='.', alpha=0.2)
polar_plot(R2, theta2, ax, lw=0.6)
ax.plot(x_tail, y_tail, '--')
ax.plot(x_tail, y_tail2, '--')

ax.set_xlim(-xmax, xmax)
ax.set_ylim(-0.2*xmax, 1.2*xmax)
ax.set_aspect('equal', adjustable='box')

phi = shell.th_infty - theta
phi1 = shell.th1_infty - th1
axx.plot(phi, phi1, alpha=0.7, label=r'$\theta_1 - \theta_{1,\infty}$ (CRW)')
axx.plot(phi, phi1, alpha=0.7, label=r'$\alpha - \theta_{1,\infty}$ (CRW)')

axx.plot(shell2.th_infty - theta2, shell2.th1_infty - th12,
         alpha=0.7, label=r'$\theta_1 - \theta_{1,\infty}$ ($k = 0$)')
axx.plot(shell2.th_infty - theta2, alph2 - shell2.th1_infty,
         alpha=0.7, label=r'$\alpha - \theta_{1,\infty}$ ($k = 0$)')

# 30 Aug 2016: plot the phi_1 = J phi + K phi^2 approximation
axx.plot(phi, fquad(phi), lw=0.5, color='k', 
         label='$J, K = {:.2f}, {:.2f}$'.format(ht.J, ht.K))
axx.plot(phi, fquad(phi, K=0.0), lw=0.5, ls='--', color='k', label=None)

axx.set_xlim(0.0, 0.8)
m = np.isfinite(phi1) & (phi < 0.8)
ymax = phi1[m].max()
print('ymax =', ymax)
axx.set_ylim(0.0, ymax)
axx.set_xlabel(r'$\theta - \theta_{\infty}$')
axx.legend(loc='upper left', fontsize='small')

axxx.plot(phi, -(phi1 - fquad(phi, K=0.0))/phi**2)
axxx.set_xlim(0.0, 0.8)
axxx.set_ylim(0.0, None)

fig.set_size_inches(4.0, 8.0)
fig.savefig(figfilename)
print(figfilename)
