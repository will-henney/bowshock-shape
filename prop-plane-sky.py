import numpy as np
import matplotlib.pyplot as plt
from equation6 import Shell

"""
In this script, we'll create the proyection of proplyd in the sky plane, assuming the plane of proplyd is rotated respect plane of sky by an angle i
"""

def OMEGA(r,rpdr,h):
    """
    Approximation to omega = (1/R) d R / d theta
    """
    dRdth = (rpdr - r)/h
    w = dR/r
    return w


def omega(r,theta):
    """
    Will's version of omega = (1/R) d R / d theta
    """
    om = np.zeros_like(r)
    om[:-1] = np.diff(np.log(r))/np.diff(theta)
    return om
    

#First, design a set of beta values to calculate R(theta) (Well, first than that, desgign a theta array)

N = 1000
h=1e-3
theta = np.linspace(0,np.pi,N)
beta = np.array([1.0, 0.2, 0.05, 0.02, 0.005])
innertype = ['isotropic']# 'proplyd']
# Use semi-open interval -pi->pi for phi
phi = np.linspace(-np.pi,np.pi,N, endpoint=False) 
#Phi is the azimutal angle in spherical coordinates to make the revolution surface of proplyd
#Meshgrid theta and phi:
THETA,PHI = np.meshgrid(theta,phi)
DPHI = 2*np.pi/N
i = np.linspace(0.0, 0.5*np.pi, 19)
for inn in innertype:
    print "******{} case******".format(inn)
    for b in beta:
        shell = Shell(beta=b,innertype=inn)
#Calculing R(theta)
        R=shell.radius(theta)
        Rdr = shell.radius(theta + h)
        R[R<=0] = np.nan # Set neagtive R to NaN
        # W = OMEGA(R,Rdr,h)
        W = omega(R, theta)

        # Check that R(theta) is working correctly
        plt.plot(theta, R)
        plt.yscale('log')
        plt.savefig("R-vs-theta-{}-{}.png".format(inn, b))
        plt.clf()

        # Check that omega is working correctly
        plt.plot(theta, W)
        plt.savefig("omega-vs-theta-{}-{}.png".format(inn, b))
        plt.clf()

        # Manage two figures at once
        xy_fig = plt.figure()     # Shape in xy plane
        xy_ax = xy_fig.add_subplot(111)
        thph_fig = plt.figure()   # Variation of ph with theta
        thph_ax = thph_fig.add_subplot(111)

        x,y,z = R*np.cos(THETA),R*np.sin(THETA)*np.cos(PHI),R*np.sin(THETA)*np.sin(PHI)
        for j in i:
            SenodePhiT = np.tan(j) * ( ( 1+ W*np.tan(THETA) )/( W-np.tan(THETA) ) )
            mask = (np.sin(PHI-0.5*DPHI) <= SenodePhiT) \
                   & (SenodePhiT <= np.sin(PHI +0.5*DPHI)) #change added
            label = "i = {:.1f}".format(np.degrees(j))
            if np.any(mask):
                # only plot a line if there really is one
                thph_ax.plot(np.degrees(THETA[mask]), np.degrees(PHI[mask]), 
                             '-', alpha=0.5, lw=3, label=label)
                xi,yi,zi = x*np.cos(j) + z*np.sin(j), y , -x*np.sin(j)+z*np.cos(j) 
                xy_ax.plot(xi[mask], yi[mask], '-', alpha=0.5, lw=3, label=label)

#first trial: plot (xi,yi) when i = 0 and phi = 0 (Proplyd in the XY plane without any rotation, we should reproduce the last cases)
#second trial: plot (xi,yi) for i = pi/6 (30 degrees)
#Third trial: plot (xi,yi) for different i values
#plt.plot(xi[PHI==0],yi[PHI==0],label = 'No rotation')
            #plt.plot(xi[PHI==np.arcsin(SenodePhiT[mask])],yi[PHI==np.arcsin(SenodePhiT[mask])],label = '30 deg rotation')
        thph_ax.legend(loc="lower right", prop=dict(size="x-small"))
        thph_ax.axis([0.0, 180.0, -100.0, 10.0])
        thph_ax.set_xlabel('theta')
        thph_ax.set_ylabel('phi')
        thph_ax.set_title('Bowshock shapes for {} wind and beta = {}'.format(inn,b))
        thph_fig.savefig('{}-wind-theta-phi-beta = {}.png'.format(inn,b))
        thph_fig.clf()

        xy_ax.legend(loc="upper right", prop=dict(size="x-small"))
        xy_ax.axis([-1.4, 0.6, -0.1, 1.9])
        xy_ax.set_xlabel('x')
        xy_ax.set_ylabel('y')
        xy_ax.set_title('Bowshock shapes for {} wind and beta = {}'.format(inn,b))
        xy_fig.savefig('{}-wind-plane-sky-beta = {}.png'.format(inn,b))
        xy_fig.clf()


#Trial succesful!!
#Next step: Try to do the same for different vales of beta and also for the proplyd case
#Now we need R90 vs R0 plots
#We nned to do everything again but with just 2 points: at theta= 0, and theta = np.pi/2, and just 1 graph

