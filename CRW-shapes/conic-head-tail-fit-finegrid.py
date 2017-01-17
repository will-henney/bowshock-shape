import sys
import numpy as np
import astropy.modeling.fitting
from astropy.table import Table
from matplotlib import pyplot as plt
import seaborn as sns
import equation6
import conic_parameters
import theta_ratio_fit 
sys.path.append('../conic-projection')
from conproj_utils import Conic

XI_LIST = [None] + list(np.linspace(1.0, 0.2, 17))
BETA_LIST = np.logspace(-4.5, -1.0, 15)
nxi, nbeta = len(XI_LIST), len(BETA_LIST)

methods = {
    'headtail': 'match head to tail',
    'gradient': 'match R90 and gradient',
    'asymptote': 'match R90 and asymptote',
}
# Method used to get initial guess at tail parameters It isn't
# critical how good this is since the fitting will sort things out
# independently of the starting point
approx_method = methods['gradient']

ntheta = 100
theta = np.linspace(0.0, np.pi, ntheta)

figfilename = sys.argv[0].replace('.py', '.pdf')

sns.set_style('whitegrid')
sns.set_color_codes('dark')

fig, axes = plt.subplots(nxi, nbeta, sharex=True, sharey=True)

xmin, xmax = -5.0, 2.1
ymin, ymax = -0.1, 7.0
# xmin, xmax = -7.0, 4.1
# ymin, ymax = -0.1, 11.0

ytop = ymin + 0.98*(ymax - ymin)
xright = xmin + 0.98*(xmax - xmin)
whitebox = {'edgecolor': 'none', 'facecolor': 'white',
            'alpha': 0.7, 'boxstyle': 'round,pad=0.1'}

outtab = {
    'beta': [],
    'xi': [],
    'x0': [],
    'a': []
}


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

        # Fit to head and approximate guess at fit to tail
        ht = conic_parameters.HeadTail(beta, xi=xi, xmin=0.0, method=approx_method)
        # Convert the tail to the form required for fitting
        model = theta_ratio_fit.hyperbola_ratio(ht.a_t, x0=ht.x0_t,
                                                tau=np.tan(ht.theta_t), D=ht.D)
        # Freeze the parameters that are external constraints
        model.tau.fixed = True
        model.D.fixed = True
        # Only fit over a certain range of angles
        #mask = (np.degrees(theta) > 100.0) & (np.degrees(theta) < 150.0)
        mask = (np.degrees(theta) > 100.0) & (x_crw > xmin)

        # Now do the fit to get an improved tail hyperbola
        best_model = fit(model, theta[mask], ratio[mask])

        # Update the tail component of the HeadTail instance
        ht.a_t = best_model.a.value
        ht.x0_t = best_model.x0.value

        # And calculate Cartesian arrays for the shapes
        x_head = ht.x_head(ht.t_h)
        y_head = ht.y_head(ht.t_h)

        x_tail = ht.x_tail(ht.t_t)
        y_tail = ht.y_tail(ht.t_t)

        # asymptote to tail
        ya = (ht.x0_t - xa)*np.tan(ht.theta_t)
        ax.plot(xa, ya, lw=0.3, color='orange')

        ax.plot(x_crw, y_crw, lw=4, color='y', alpha=0.7)
        ax.plot(x_head, y_head, '--', color='g')
        ax.plot(x_tail, y_tail, '-', dashes=[8, 4, 2, 4], color='r')
        ax.plot(0.0, 0.0, 'o', color='k')
        ax.plot(ht.x0_t - ht.a_t, 0.0, '.', color='k')
        ax.axhline(ls=':')
        ax.axvline(ls=':')
        ax.set_aspect('equal', adjustable='box-forced')

        if xi is None:
            text = r'Isotropic'
        else:
            text = r'Anisotropic, $\xi = {:.1f}$'.format(xi)
        text += '\n' + r'$\beta = {:.4f}$'.format(beta)
        text += ', ' + r'$D = {:.1f}$'.format(ht.D)
        text += '\n' + r'$x_t = {:.1f}$'.format(ht.x0_t)
        text += ', ' + r'$x_t - a_t = {:.1f}$'.format(ht.x0_t - ht.a_t)
        text += ', ' + r"$\theta_t = {:.1f}$".format(np.degrees(-ht.theta_t))
        text += '\n' + r'$x_h = {:.1f}$'.format(ht.x0_h)
        text += ', ' + r"$\theta_h = {:.1f}$".format(np.degrees(ht.theta_h))
        ax.text(xright, ytop, text,
                ha='right', va='top', bbox=whitebox, fontsize='small')

        print('beta = {:.4f}, xi = {:.2f}, x0 = {:.2f}, a = {:.2f}'.format(
            beta, xi or 0.0, ht.x0_t, ht.a_t))
        outtab['beta'].append(beta)
        outtab['xi'].append(xi)
        outtab['x0'].append(ht.x0_t)
        outtab['a'].append(ht.a_t)

# Put axis labels on lower left panel only
axes[-1, 0].set(
    xlim=[xmin, xmax], ylim=[ymin, ymax],
    xlabel=r'$x / r_{0}$', ylabel=r'$y / r_{0}$',
)

fig.set_size_inches(3*nbeta, 3*nxi)
fig.tight_layout()
fig.savefig(figfilename)
print(figfilename)

tabfilename = figfilename.replace('.pdf', '.tab')
Table(outtab).write(tabfilename, format='ascii.tab')
