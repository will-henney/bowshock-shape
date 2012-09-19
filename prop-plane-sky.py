import numpy as np
import matplotlib.pyplot as plt
from equation6 import Shell

"""
In this script, we'll create the proyection of proplyd in the sky plane, assuming the plane of proplyd is rotated respect plane of sky by an angle i
"""

def OMEGA(r,rpdr,h):

    dR = (rpdr - r)/h
    w = dR/r
    return w

#First, design a set of beta values to calculate R(theta) (Well, first than that, desgign a theta array)

N = 50
h=1e-6
theta = np.linspace(0,np.pi,N)
beta = np.array([1.0, 0.2, 0.05, 0.02, 0.005])
innertype = ['isotropic','proplyd']
phi = np.linspace(0,2*np.pi,N) #Phi is the azimutal angle in spherical coordinates to make the revolution surface of proplyd
#Meshgrid theta and phi:
THETA,PHI = np.meshgrid(theta,phi)
DPHI = 2*np.pi/N
i = np.array([0.0,np.pi/6,np.pi/4,np.pi/3])
for inn in innertype:
    print "******{} case******".format(inn)
    for b in beta:
        shell = Shell(beta=b,innertype=inn)
#Calculing R(theta)
        R=shell.radius(theta)
        Rdr = shell.radius(theta + h)
        R[R<=0] = np.nan # Set neagtive R to NaN
        W = OMEGA(R,Rdr,h)
        x,y,z = R*np.cos(THETA),R*np.sin(THETA)*np.cos(PHI),R*np.sin(THETA)*np.sin(PHI)
        for j in i:
            SenodePhiT = np.tan(j) * ( ( 1+ W*np.tan(THETA) )/( W-np.tan(THETA) ) )
            mask = np.abs(SenodePhiT) <= 1.  #Discard non real values for the angle of tangent line
            mask2 = (np.sin(PHI-0.5*DPHI) <= SenodePhiT)  &  (SenodePhiT <= np.sin(PHI +0.5*DPHI)) #change added
            xi,yi,zi = x*np.cos(j) + z*np.sin(j), y , -x*np.sin(j)+z*np.cos(j) 
            plt.plot(xi[mask & mask2],yi[mask & mask2],'+',label = "i = {} (radians)".format(j))
#first trial: plot (xi,yi) when i = 0 and phi = 0 (Proplyd in the XY plane without any rotation, we should reproduce the last cases)
#second trial: plot (xi,yi) for i = pi/6 (30 degrees)
#Third trial: plot (xi,yi) for different i values
#plt.plot(xi[PHI==0],yi[PHI==0],label = 'No rotation')
#plt.plot(xi[PHI==np.arcsin(SenodePhiT[mask])],yi[PHI==np.arcsin(SenodePhiT[mask])],label = '30 deg rotation')
        plt.legend()
        plt.axis([-0.6, 0.6, -0.1, 1.1])
        plt.xlabel('z')
        plt.ylabel('r')
#plt.title('No rotation')
#plt.title('30 deg rotation')
        plt.title('Bowshock shapes for {} wind and beta = {}'.format(inn,b))
#plt.savefig('trial.png')
#plt.savefig('trial2.png')
        plt.savefig('{}-wind-plane-sky-beta = {}.png'.format(inn,b))
        plt.clf()
#Trial succesful!!
#Next step: Try to do the same for different vales of beta and also for the proplyd case
#Now we need R90 vs R0 plots
#We nned to do everything again but with just 2 points: at theta= 0, and theta = np.pi/2, and just 1 graph

