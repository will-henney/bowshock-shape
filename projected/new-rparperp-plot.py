"""
New program to plot R_par and R_perp, but divided by D

This is so that we can compare them with Jorge's plots and the figure in the Robberto paper

The tables that are written out by rparperp.py are normalized so that Rpar = 1 for i = 0
"""
from os.path import splitext, basename
import numpy as np
from matplotlib import pyplot as plt
from bowfuncs import theta_lim, theta_par_approx, \
    Rpar, Rperp, theta_par_approx, theta_par, theta_perp

betalist = np.logspace(1.e-3, 0.3, 100)

for beta in betalist: 
    # Calcluate the radius of stagnation point R0 in units of separation D
    R0_over_D = np.sqrt(beta)/(1.0 + np.sqrt(beta)) # CRW Eq (27)

    betas, incs, thpars, thperps, rpars, rperps = \
	np.loadtxt("rparperp-B%.2e.dat" % (beta), skiprows=1, unpack=True)

    rpars *= R0_over_D
    rperps *= R0_over_D
    incs  = np.degrees(incs)

    plt.plot(rpars, rperps, "o", label="beta = {:.1e}".format(beta))

plt.ylim(0.0, 4.0)
plt.xlabel("R_par")
plt.ylabel("R_perp")
plt.legend(loc="upper left")

# make sure that the PDF file has a name that matches this python script
prefix, _ = splitext(basename(__file__))
plt.savefig("{}.pdf".format(prefix))
    
