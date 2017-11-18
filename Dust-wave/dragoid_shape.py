import sys
import numpy as np
from astropy.table import Table
import statsmodels.api as sm
from bow_projection import Spline_R_theta_from_grid

class Dragoid(object):
    def __init__(self, alpha, mu=None, lowess_frac=None):
        if mu is None:
            astring = 'dust-couple-stream'
        else:
            astring = 'dust-couple-div-stream'
        astring += f'-alpha{int(100*alpha):03d}'
        if mu is not None:
            astring += f'-mu{int(100*mu):03d}'
        astring += '.tab'
        self.label = fr"$\alpha_\mathrm{{drag}} =  {alpha:.02f}$"
        if mu is not None:
            self.label += ', ' + fr"$\mu =  {mu:.02f}$"
        t = Table.read(astring, format='ascii.tab')
        dth = np.pi/len(t)
        self.thgrid = t['theta'] + 0.5*dth
        self.Rgrid = t['R']/t['R'][0]
        self.thgrid = np.concatenate([-self.thgrid[::-1], self.thgrid])
        self.Rgrid = np.concatenate([self.Rgrid[::-1], self.Rgrid])
        if lowess_frac is not None:
            # Optionally smooth the shape before fitting spline
            self.Rgrid = sm.nonparametric.lowess(
                self.Rgrid, self.thgrid, frac=lowess_frac,
                is_sorted=True, return_sorted=False)
        self.splinefit = Spline_R_theta_from_grid(
              theta_grid=self.thgrid, R_grid=self.Rgrid)

    def __call__(self, theta):
        # When called as a function, give the spline fitted result
        return self.splinefit(theta)

if __name__ == "__main__":

    from matplotlib import pyplot as plt
    import seaborn as sns

    lib_name = sys.argv[0].replace('.py', '')
    figfile = f"test_{lib_name}_radius.pdf"

    sns.set_style('ticks')
    fig, ax = plt.subplots()

    th = np.linspace(-np.pi, np.pi, 1001)
    th_dg = np.degrees(th)

    alphas = [0.25, 0.5, 1.0, 2.0] + [4.0, 4.0]
    mus = [None]*4 + [0.2, 0.8]
    for alpha, mu in zip(alphas, mus):
        shape = Dragoid(alpha=alpha, mu=mu, lowess_frac=None)
        ax.plot(np.degrees(shape.thgrid), shape.Rgrid,
                color='b', alpha=0.2, lw=2, label='_nolabel_')
        ax.plot(th_dg, shape(th), lw=0.8, label=shape.label)

    ax.legend(title=r"Dragoid shapes")
    ax.set(
        xlabel=r"Polar angle: $\theta$, degrees",
        ylabel=r"$R$",
        xlim=[0, 180],
        yscale='log',
        ylim=[0.9, 200.0],
        xticks=[0, 30, 60, 90, 120, 150, 180],
    )
    sns.despine()
    fig.tight_layout()
    fig.savefig(figfile)
    print(figfile, end='')
