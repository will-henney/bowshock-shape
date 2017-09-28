import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker
import seaborn as sns
import conic_parameters

try: 
    xiset = sys.argv[1]
    plotfile = sys.argv[0].replace('.py', f'-{xiset}.pdf')
    assert xiset in 'ab'
    istart = -2 if xiset == 'a' else -1
except:
    sys.exit(f"Usage: {sys.argv[0]} a|b")

sns.set_style('white')
fig, ax = plt.subplots(figsize=(5, 5))


left_annotate_pars = dict(xytext=(-5, 5), ha='right', va='bottom')
right_annotate_pars = dict(xytext=(5, -5), ha='left', va='top')


Rc_grid = np.linspace(0.0, 10.0, 2000)
R90_T0_grid = np.sqrt(2*Rc_grid)
R90_T1_grid = np.sqrt(2*Rc_grid - 1.0)
R90_T1_grid[~np.isfinite(R90_T1_grid)] = 0.0 

ax.fill_between(Rc_grid, R90_T1_grid, R90_T0_grid, color='k', alpha=0.2)
ax.fill_between(Rc_grid, R90_T0_grid, color='k', alpha=0.1)
ax.plot(Rc_grid, R90_T0_grid, c='k', lw=0.5)
ax.axhline(1.0, lw=0.5, alpha=0.5, color='k', zorder=-1)
ax.axvline(1.0, lw=0.5, alpha=0.5, color='k', zorder=-1)
ax.plot([0.0, 10.0], [0.0, 10.0], lw=0.5, alpha=0.5, color='k', zorder=-1)

inc = np.linspace(0.0, 0.5*np.pi, 500, endpoint=False)
inc_deg = np.degrees(inc)


XI_LIST = [None, 1.0, 0.8, 0.4]
BETA_LIST = [0.2, 0.1, 0.05, 0.02, 0.005, 1e-6]
nxi, nbeta = len(XI_LIST), len(BETA_LIST)

cols = sns.color_palette('magma', n_colors=nxi)
annot_pars_list = [left_annotate_pars]*2 + [right_annotate_pars]*2 
for beta in BETA_LIST[::-1]:
    for xi, col, annot_pars in list(zip(XI_LIST, cols, annot_pars_list))[istart::-2]:

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

        # if Rc == 1.0:
        #     label = fr'$\xi = {xi:.1f}$'
        # else:
        #     label = None
        if beta == BETA_LIST[0]:
            label = 'CRW' if xi is None else fr'$\xi = {xi:.1f}$'
        else:
            label = None

        # Find minimum difference between head and tail values of R90
        dR = np.abs(R90_h_prime - R90_t_prime)
        mm = np.isfinite(dR)
        i0 = np.argmin(dR[mm])


        # Masks for high and low inclinations (overlapping range)
        overlap_deg = 0.0
        mlo = inc_deg <= inc_deg[i0] + overlap_deg
        mhi = inc_deg >= inc_deg[i0] - overlap_deg
        # mlo = inc_deg <= 40.0
        # mhi = inc_deg >= 30.0

        # Ensure that LOS is not inside the tail cone
        mhi = mhi & (inc < 0.5*np.pi - ht.theta_t)
        # Plot the head discriminant for low inclinations
        ax.plot(tilde_Rc_h_prime[mlo], R90_h_prime[mlo]/R0_h_prime[mlo],
                c=col, label=None, lw=1.6, alpha=0.7)
        # Put a dot at the i=0 case
        ax.plot([tilde_Rc_h], [R90_h/R0_h], '.', c=col, alpha=1.0)
        # Plot the combined discriminant for high inclinations
        ax.plot(tilde_Rc_h_prime[mhi], R90_t_prime[mhi]/R0_h_prime[mhi],
                c=col, label=label, lw=0.8, alpha=1.0)
        # Label the dot with the cross-over inclination
        ax.annotate(rf'$\beta = \mathrm{{{beta:g}}}$, $i^*\!\!\! = {inc_deg[i0]:.0f}^\circ$',
                    xy=(tilde_Rc_h, R90_h/R0_h),
                    textcoords='offset points',
                    fontsize='x-small', color=col, **annot_pars)
        # Put a cross at the Wilkinoid coordinates: [5/3, sqrt(3)]
        ax.plot([5./3.], [np.sqrt(3.0)], '+', c='w', ms=10, alpha=1.0)

ax.legend(ncol=1, fontsize='small', frameon=True)
ax.set(
    yscale='linear',
    xscale='linear',
    xlim=[0.0, 5.1],
    ylim=[0.0, 5.1],
#    ylim=[-3.0, 1.1],
    xlabel=r"Projected dimensionless radius of curvature: $\widetilde{R}_{c}{}'$",
    ylabel=r"Projected dimensionless perpendicular radius: $\widetilde{R}_{90}{}'$",
)        


fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
