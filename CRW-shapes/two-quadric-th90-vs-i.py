import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker
import seaborn as sns
import conic_parameters

plotfile = sys.argv[0].replace('.py', '.pdf')

sns.set_style('white')
fig, (axx, ax) = plt.subplots(2, 1, sharex=True, figsize=(5, 5))



inc = np.linspace(0.0, 0.5*np.pi, 500, endpoint=False)
inc_deg = np.degrees(inc)


XI_LIST = [None, 1.0, 0.8, 0.4]
BETA_LIST = [0.2, 0.1, 0.05, 0.02, 0.005, 1e-6]
nxi, nbeta = len(XI_LIST), len(BETA_LIST)

dashes_solid = []
dashes_dashed = [3, 2]
dashes_dotted = [1, 2]
dashes_dot_dashed = [1, 2, 4, 2]
dashes_triple_dot_dashed = [1, 2, 1, 2, 1, 2, 4, 2]
dashes_triple_dot_spaced = [1, 2, 1, 2, 1, 6]
dashes_styles = [dashes_solid, dashes_dashed, dashes_dotted,
                 dashes_dot_dashed, dashes_triple_dot_dashed,
                 dashes_triple_dot_spaced]

cols = sns.color_palette('magma', n_colors=nxi)

ax.axhspan(100.0, 110.0, facecolor='k', alpha=0.1)

for beta, dashes in list(zip(BETA_LIST, dashes_styles))[::-1]:
    for xi, col in list(zip(XI_LIST, cols))[-4::-1]:

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

        t2i = np.tan(inc)**2
        tan_th90 = -np.sqrt(t2i*(2.0 + T_h*t2i) + (2.0 - T_h/tilde_Rc_h)/tilde_Rc_h)/t2i
        th90_h = 180.0 + np.degrees(np.arctan(tan_th90))


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

        tan_th90 = -np.sqrt(t2i*(2.0 + T_t*t2i) + (2.0 - T_t/tilde_Rc_t)/tilde_Rc_t)/t2i
        th90_t = 180.0 + np.degrees(np.arctan(tan_th90))

        # if Rc == 1.0:
        #     label = fr'$\xi = {xi:.1f}$'
        # else:
        #     label = None
        if beta == BETA_LIST[0]:
            xilabel = 'CRW' if xi is None else fr'$\xi = {xi:.1f}$'
        else:
            label = None

        # Find minimum difference between head and tail values of R90
        dR = np.abs(R90_h_prime - R90_t_prime)
        mm = np.isfinite(dR)
        i0 = np.argmin(dR[mm])


        # Masks for high and low inclinations (overlapping range)
        mlo = inc_deg <= inc_deg[i0]
        mhi = inc_deg >= inc_deg[i0]
        # Ensure that LOS is not inside the tail cone
        mhi = mhi & (inc < 0.5*np.pi - ht.theta_t)


        ax.plot(inc_deg[mlo], th90_h[mlo],
                c=col, label=None, lw=1.5, dashes=dashes, alpha=0.8)
        ax.plot(inc_deg[mhi], th90_t[mhi],
                c=col, label=label, lw=0.7, dashes=dashes, alpha=0.8)

        axx.plot(inc_deg[mlo], R90_h_prime[mlo],
                 c=col, label=None, lw=1.5, dashes=dashes, alpha=0.8)
        axx.plot(inc_deg[mhi], R90_t_prime[mhi],
                 c=col, label=label, lw=0.7, dashes=dashes, alpha=0.8)
        axx.plot([inc_deg[i0]], [R90_t_prime[i0]], '.', c=col, label=None) 

ax.legend(ncol=1, fontsize='xx-small', frameon=True)
ax.set(
    yscale='linear',
    xscale='linear',
    xlim=[0.0, 90.0],
    ylim=[80.0, 180.0],
#    ylim=[-3.0, 1.1],
    xlabel=r'Inclination, degrees',
    ylabel=r"Body-frame polar angle of perpendicular projected axis: $\theta_{90}$, degrees",
    xticks=[15, 30, 45, 60, 75, 90],
)        
axx.set(
    yscale='log',
    xscale='linear',
    ylim=[1.0, 100.0],
    ylabel=r"Projected perpendicular radius, $R'_{90}$",
)        

sns.despine()
#fig.tight_layout()
fig.savefig(plotfile)
print(plotfile, end='')
