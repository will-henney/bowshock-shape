import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker
import seaborn as sns
import conic_parameters

plotfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('white')
fig, axes = plt.subplots(3, 3, figsize=(9, 9), sharex=True, sharey=True)

incs_deg = 10.0*np.arange(9)

nbeta = 30
#betas = np.logspace(-5.0, -0.5, nbeta)
betas = np.linspace(0.003, 0.5, nbeta)**2
XI_LIST = [None, 1.0, 0.8, 0.4]
nxi = len(XI_LIST)

Rc_grid = np.linspace(0.0, 10.0, 2000)
R90_T0_grid = np.sqrt(2*Rc_grid)
R90_T1_grid = np.sqrt(2*Rc_grid - 1.0)
R90_T1_grid[~np.isfinite(R90_T1_grid)] = 0.0 


cols = sns.color_palette('magma', n_colors=nxi)
for ax, inc_deg in zip(axes.flat, incs_deg):
    ax.fill_between(Rc_grid, R90_T1_grid, R90_T0_grid, color='k', alpha=0.2)
    ax.fill_between(Rc_grid, R90_T0_grid, color='k', alpha=0.1)
    ax.plot(Rc_grid, R90_T0_grid, c='k', lw=0.5)
    ax.axhline(1.0, lw=0.5, alpha=0.5, color='k', zorder=-1)
    ax.axvline(1.0, lw=0.5, alpha=0.5, color='k', zorder=-1)
    ax.plot([0.0, 10.0], [0.0, 10.0], lw=0.5, alpha=0.5, color='k', zorder=-1)
    for xi, col in list(zip(XI_LIST, cols)):
        for beta in betas:
            # Fit to head and analytic fit to fit to tail
            ht = conic_parameters.HeadTail(beta, xi=xi, xmin=0.0, method='analytic fit')
            # Parameters for head conic
            T_h = ht.sig_h*ht.tau_h**2
            tilde_Rc_h = ht.A_h
            R0_h = 1.0
            R90_h = ht.R90

            # Parameters for tail conic
            T_t = -ht.tau_t**2
            R0_t = ht.x0_t - ht.a_t
            # Equation E from notes
            tilde_Rc_t = np.abs(T_t)*ht.a_t/R0_t
            R90_t = R0_t * np.sqrt(2*tilde_Rc_t - T_t)
            T_combine = 2*tilde_Rc_h - (R90_t / R0_h)**2

            inc = np.radians(inc_deg)

            # Projected head quantities as functions of inc
            f_h = np.sqrt(1.0 + T_h * np.tan(inc)**2)
            tilde_Rc_h_prime = tilde_Rc_h / (
                np.cos(inc)**2 * f_h * (
                    1.0 + (tilde_Rc_h / T_h) * (f_h - 1.0) 
                )
            )
            T_h_prime = T_h / (np.cos(inc)**2 * f_h**2)
            R0_h_prime = R0_h * np.cos(inc) * (
                1.0 + (tilde_Rc_h / T_h) * (f_h - 1.0)
            )
            R90_h_prime = R0_h_prime * np.sqrt(2*tilde_Rc_h_prime - T_h_prime)


            # Projected tail quantities as functions of inc
            f_t = np.sqrt(1.0 + T_t * np.tan(inc)**2)
            # Equation B from notes
            T_t_prime = T_t / f_t**2 / np.cos(inc)**2
            # Equation D from notes
            R0_t_prime = R0_t * np.cos(inc) * (
                1.0 + (tilde_Rc_t / T_t) * (f_t - 1.0)
            )
            # Equation C from notes
            tilde_Rc_t_prime = tilde_Rc_t / (
                np.cos(inc)**2 * f_t * (
                    1.0 + (tilde_Rc_t / T_t) * (f_t - 1.0) 
                )
            )
            # Equation A from notes
            R90_t_prime = R0_t_prime * np.sqrt(2*tilde_Rc_t_prime - T_t_prime)

            # Finally, the combined discriminant (equation F from notes)
            T_combine_prime = 2*tilde_Rc_h_prime - (R90_t_prime / R0_h_prime)**2

            if inc_deg < 30.0:
                # Plot the head for low inclinations
                y = R90_h_prime/R0_h_prime
            else:
                # Plot the tail for high inclinations
                y = R90_t_prime/R0_h_prime
            ax.scatter([tilde_Rc_h_prime], [y],
		       c=col, edgecolors='none',
		       marker='o', s=25*R0_h_prime/R0_h, alpha=0.4)

            ax.text(3.0, 0.5, rf'$i = {inc_deg:.0f}^\circ$',
                    bbox={'facecolor': 'w', 'alpha': 0.8, 'edgecolor': 'none'})


axes[-1, 0].set(
    yscale='linear',
    xscale='linear',
    xlim=[0.0, 5.1],
    ylim=[0.0, 5.1],
    xlabel=r"$\widetilde{R}_{c}{}'$",
    ylabel=r"$\widetilde{R}_{90}{}'$",
)        


fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
