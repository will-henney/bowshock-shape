import sys
import numpy as np
import astropy.modeling.fitting
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.interpolate import interp1d
import equation6
import conic_parameters
import theta_ratio_fit 
sys.path.append('../conic-projection')
from conproj_utils import Conic

XI_LIST = [None, 1.0, 0.8, 0.4]
BETA_LIST = [0.1, 0.01, 0.0001]
nxi, nbeta = len(XI_LIST), len(BETA_LIST)

ntheta = 100
theta = np.linspace(0.0, np.pi, ntheta)

figfilename = sys.argv[0].replace('.py', '.pdf')

sns.set_style('whitegrid')
sns.set_color_codes('dark')

NROWS = 2
fig, axes = plt.subplots(nxi, nbeta, sharex=True, sharey=True)

xmin, xmax = -5.0, 2.1
ymin, ymax = -0.1, 7.0
# xmin, xmax = -7.0, 4.1
# ymin, ymax = -0.1, 11.0

ytop = ymin + 0.98*(ymax - ymin)
xright = xmin + 0.98*(xmax - xmin)
whitebox = {'edgecolor': 'none', 'facecolor': 'white',
            'alpha': 0.7, 'boxstyle': 'round,pad=0.1'}


# x-data for tail asymptote
xa = np.linspace(xmin, xmax, 2)

# Set up fitter for fitting the tail
fit = astropy.modeling.fitting.LevMarLSQFitter()

for j, xi in enumerate(XI_LIST):
    for i, beta in enumerate(BETA_LIST[::-1]):
        ax = axes[j, i]

        # The exact solution to the shell
        if xi is None:
            shell = equation6.Shell(innertype='isotropic', beta=beta)
        else:
            shell = equation6.Shell(innertype='anisotropic', beta=beta, xi=xi)
        R, theta1 = shell.radius(theta, full=True)
        ratio = theta1/theta
        R_crw = R/shell.R0
        x_crw = R_crw*np.cos(theta)
        y_crw = R_crw*np.sin(theta)

        # Fit to head and analytic fit to fit to tail
        ht = conic_parameters.HeadTail(beta, xi=xi,
				       xmin=0.0, method='analytic fit')

        # And calculate Cartesian arrays for the shapes
        x_head = ht.x_head(ht.t_h)
        y_head = ht.y_head(ht.t_h)
        x_tail = ht.x_tail(ht.t_t)
        y_tail = ht.y_tail(ht.t_t)

        # Work in the ratios theta_1 / theta ...
        # ... first for the head ...
        theta_head = np.arctan2(y_head, x_head)
        theta1_head = np.arctan2(y_head, ht.D - x_head)
        ratio_head_func = interp1d(theta_head, theta1_head/theta_head,
                                   fill_value='extrapolate')
        ratio_head = ratio_head_func(theta)
        # ... and then for the tail
        model = theta_ratio_fit.hyperbola_ratio(ht.a_t, x0=ht.x0_t,
                                                tau=np.tan(ht.theta_t), 
                                                D=ht.D)


        # Fractional residual in theta_1 is same as fractional
        # residual in the ratio
        resid_head = (ratio_head - ratio)/ratio
        resid_tail = (model(theta) - ratio)/ratio

        # Heuristic to find the switc-over point between the head and the tail
        m = (theta < 0.2) | (resid_head**2 < resid_tail**2)
        thm = np.min(theta[~m])
        m[theta > thm] = False

        ax.plot(np.degrees(theta[~m][:-1]), resid_tail[~m][:-1],
                c='r', label='Tail fit')
        ax.plot(np.degrees(theta[:-1]), resid_tail[:-1],
                c='r', lw=1, ls=':', label=None)
        ax.plot(np.degrees(theta[m]), resid_head[m],
                c='g', label='Head fit')
        ax.plot(np.degrees(theta), resid_head,
                c='g', lw=1, ls=':', label=None)
        ax.axvline(np.degrees(shell.th_infty), lw=0.5,
                   c='k', ls='--')

        if xi is None:
            text = r'Isotropic'
        else:
            text = r'Anisotropic, $\xi = {:.1f}$'.format(xi)
        text += '\n' + r'$\beta = {:.4f}$'.format(beta)
        ax.text(0.02, 0.98, text,
                ha='left', va='top', transform=ax.transAxes,
                bbox=whitebox, fontsize='small')


# Put legend on upper left panel only
axes[0, 0].legend(fontsize='x-small', frameon=True)
# Put axis labels on lower left panel only
axes[-1, 0].set(
    xlim=[0.0, 180.0], ylim=[-0.25, 0.25],
    xticks=[0, 30, 60, 90, 120, 150, 180],
    yticks=[-0.2, -0.1, 0.0, 0.1, 0.2],
    xlabel=r'$\theta$, degrees', ylabel=r'$\Delta \theta_1 / \theta_{1}$', 
)

fig.set_size_inches(2*nbeta, 1.5*nxi)
fig.tight_layout()
fig.savefig(figfilename)
print(figfilename, end='')
