"Calculate parallel and perpendicular projected bowshock radii"
from bowfuncs import *

import numpy as np

# betalist = [0.3, 0.03, 0.003]
# betalist = [0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 3.e-4, 1.e-4, 3.e-5, 1.e-5]
betalist = [0.3, 0.1, 0.01, 0.001, 1.e-4, 1.e-5, 1.e-7]

n = 200
deg = 180.0/np.pi
ratiolines = []
rparlines = []
rperplines = []
thparlines = []
thparapproxlines = []
verbose = True
for beta in betalist:
    f = open("rparperp-new-B%.2e.dat" % (beta), "w")
    f.write(
	"%s\t"*6 % ("beta", "i", "thpar", "thperp", "rpar", "rperp") + "\n"
	)
    thlim = theta_lim(beta)
    incs = np.linspace(0.0, 0.98*(thlim - 0.5*pi), n)
    for i, inc in enumerate(incs):
        rpar, rperp, thpar, thperp = find_Rpar_Rperp(beta, inc, full=True)
	f.write(
	    "%.3e\t"*6 % (beta, inc, thpar, thperp, rpar, rperp) + "\n"
	    )
    f.close()

