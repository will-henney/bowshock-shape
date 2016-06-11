import numpy as np
import equation6
from astropy.table import Table, Column

def find_alpha(beta, xi, dtheta=0.01, ntheta=10, npoly=2):
    """Find the coefficient of theta**2 in the Taylor expansion of R(theta)"""
    theta = dtheta*np.arange(ntheta)
    thsq = theta**2
    shell = equation6.Shell(innertype='anisotropic', beta=beta, xi=xi)
    R = shell.radius(theta)
    R0 = np.sqrt(beta)/(1 + np.sqrt(beta))
    ahat = (R/R0 - 1) / thsq
    c = np.polyfit(thsq[2:], ahat[2:], npoly)
    return c[-1]



if __name__ == '__main__':

    xi_range = np.linspace(0.1, 1.0, 10)
    beta_range = np.linspace(0.001, 1.0, 300)

    tab = Table()
    tab.add_column(Column(name='beta', data=beta_range, format='{:.6f}'))
    for xi in xi_range:
        # print(xi)
        tab.add_column(Column(name='xi = {:.1f}'.format(xi),
                              data=[find_alpha(beta, xi) for beta in beta_range],
                              format='{:.6f}'))

    tabfilename = 'alpha-vs-beta-xi.tab'
    tab.write(tabfilename, format='ascii.tab')
    print(tabfilename)
