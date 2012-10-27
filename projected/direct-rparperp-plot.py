"""
New program to plot R_par and R_perp, but divided by D

This is so that we can compare them with Jorge's plots and the figure in the Robberto paper

The tables that are written out by rparperp.py are normalized so that Rpar = 1 for i = 0
"""
from os.path import splitext, basename
import numpy as np
from matplotlib import pyplot as plt
from bowfuncs import theta_lim, find_Rpar_Rperp, R0


tol = 0.01

# First plot solutions as a function of inc for constant beta
# betas = np.logspace(-4.0, -0.2, 10)
# for beta in betas: 
#     thlim = theta_lim(beta)
#     print "beta = ", beta, "thlim = ", np.degrees(thlim)
#     inc_lim = (1.0 - tol)*(thlim - 0.5*np.pi)
#     incs = np.linspace(0.0, inc_lim, 101)
#     rpars, rperps = list(), list()
#     for inc in incs:
#         rpar, rperp = find_Rpar_Rperp(beta, inc)
#         rpars.append(R0(beta)*rpar)
#         rperps.append(R0(beta)*rperp)

#     plt.plot(rpars, rperps, "--", label="beta = {:.4f}".format(beta))

# Second, plot solutions as function of beta for constant inc
incs_deg = np.linspace(0.0, 90.0, 18, endpoint=False)
betas = np.logspace(-6.0, -0.2, 300)
for inc_deg in incs_deg:
    inc = np.radians(inc_deg)
    print "inc = ", inc_deg
    rpars, rperps = list(), list()
    for beta in betas:
        thlim = theta_lim(beta)
        inc_lim = (1.0 - tol)*(thlim - 0.5*np.pi)
        if inc > inc_lim:
            break
        rpar, rperp = find_Rpar_Rperp(beta, inc)
        rpars.append(R0(beta)*rpar)
        rperps.append(R0(beta)*rperp)

    plt.plot(rpars, rperps, "-", label="inc = {:.0f}".format(inc_deg))
    

plt.xlim(0.0, 0.35)
plt.ylim(0.0, 1.0)
plt.xlabel("R_par")
plt.ylabel("R_perp")
plt.legend(loc="upper left", ncol=2, prop=dict(size="xx-small"))

# make sure that the PDF file has a name that matches this python script
prefix, _ = splitext(basename(__file__))
plt.savefig("{}.pdf".format(prefix))
    
