import numpy as np
from scipy.optimize import brentq
from scipy.misc import derivative

# This should be around sqrt(eps), where eps is machine precision
DX_FOR_NUMERICAL_DERIVATIVE = 3.0*np.sqrt(np.finfo(1.0).resolution)

def omega(func_R, theta, *args_for_func_R):
    """Find omega = R^{-1} d R/d theta 

Note that theta may be an array. Any extra arguments are passed to
`func_R` after `theta`

    """
    def log_R(theta, *args):
        return np.log(func_R(theta, *args))

    return derivative(log_R, theta,
                      dx=DX_FOR_NUMERICAL_DERIVATIVE, args=args_for_func_R)



# * Example shapes 

def wilkinoid_R_theta(theta):
    """Wilkin solution for wind-stream interaction"""
    return np.sqrt(3*(1.0 - theta/np.tan(theta)))/np.sin(theta)


def cantoid_R_theta(theta, beta):
    """Cantoid solution from CRW for wind-wind interaction

Returns R(theta), normalized to the stagnation radius. Extra parameter
`beta` is the momentum ratio of the two winds

    """

    # Approximate solution for theta_1, the polar angle measured from
    # the "other" star
    theta1 = np.sqrt(7.5*(-1.0 + np.sqrt(
        1.0 + 0.8*beta*(1.0 - theta/np.tan(theta)))))

    # On-axis (theta = 0) radius to stagnation point, in units of
    # star-star separation D
    R0 = np.sqrt(beta)/(1.0+np.sqrt(beta))

    # Return radius in units of R0
    return  np.sin(theta1) / np.sin(theta + theta1) / R0


def paraboloid_R_theta(theta):
    """This is the special parabola with Rc = 2"""
    return 2.0 / (1.0 + np.cos(theta))


def paraboloid_omega_true(theta):
    """Analytic omega for special parabola"""
    return np.sin(theta)  / (1.0 + np.cos(theta))


if __name__ == "__main__":
    import sys
    from matplotlib import pyplot as plt
    import seaborn as sns

    lib_name = sys.argv[0].replace('.py', '')

    sns.set_style('white')
    fig, ax = plt.subplots()

    th = np.linspace(0.0, np.pi, 501)
    th_dg = np.degrees(th)
    ax.plot(th_dg, omega(paraboloid_R_theta, th),
            label="paraboloid")
    ax.plot(th_dg, omega(wilkinoid_R_theta, th),
            label="wilkinoid")
    for beta in 0.001, 0.01, 0.1:
        ax.plot(th_dg, omega(cantoid_R_theta, th, beta),
                label=fr"cantoid $\beta = {beta:.3f}$")
    ax.legend()
    ax.axhline(1.0, xmin=0.35, xmax=0.65, color='white', lw=4, zorder=100)
    ax.axhline(1.0, xmin=0.35, xmax=0.65, color='k', lw=1, ls=':', zorder=101)
    ax.set_yscale('symlog', linthreshy=1.0, linscaley=0.5)
    ax.set(
        xlabel=r"Polar angle: $\theta$, degrees",
        ylabel=r"$\omega = R^{-1} d R / d \theta$",
        xlim=[0, 180],
        ylim=[0.0, 1e2],
        xticks=[0, 30, 60, 90, 120, 150, 180],
    )
    sns.despine()
    fig.tight_layout()
    figfile = f"test_{lib_name}_omega.pdf"
    fig.savefig(figfile)
    print(figfile, end='')
