import numpy as np
import matplotlib.pyplot as plt
from equation6 import Shell
import argparse
"""
In this script, we'll create the proyection of proplyd in the sky plane, assuming the plane of proplyd is rotated respect plane of sky by an angle i
"""

def OMEGA(r,rpdr,h):

    dR = (rpdr - r)/h
    w = dR/r
    return w

def omega(r,t):
    """
    Will's version of omega(theta), looks like works better
    uses d(log(r))/dtheta = 1/r * dr/dtheta (ingenious)
    """
    womega = np.zeros_like(r)
    womega[:-1] = np.diff(np.log(r))/np.diff(theta)
    return womega

#First, design a set of beta values to calculate R(theta) (Well, first than that, desgign a theta array)
parser = argparse.ArgumentParser(description="Rotation angle")
parser.add_argument("--iangle", "-i", default=0.0, type=float, help="Rotation of z axis in radians")
cmdargs = parser.parse_args()
i = cmdargs.iangle

N = 50
h=1e-6
theta = np.linspace(0,np.pi,N)
beta = np.array([0.2, 0.05, 0.02, 0.005])
innertype = ['isotropic','proplyd']
phi = np.linspace(-np.pi,np.pi,N,endpoint = False) #Phi is the azimutal angle in spherical coordinates to make the revolution surface of proplyd
#Meshgrid theta and phi:
THETA,PHI = np.meshgrid(theta,phi)
DPHI = 2*np.pi/N
#i = np.array([0.0,np.pi/6,np.pi/4,np.pi/3])
for inn in innertype:
    print "******{} case******".format(inn)
    for b in beta:
        shell = Shell(beta=b,innertype=inn)
#Calculing R(theta)
        R=shell.radius(theta)
#        Rdr = shell.radius(theta + h)
        R[R<=0] = np.nan # Set neagtive R to NaN
#        W = OMEGA(R,Rdr,h)
        w = omega(R,theta)
        x,y,z = R*np.cos(THETA),R*np.sin(THETA)*np.cos(PHI),R*np.sin(THETA)*np.sin(PHI)
        SenodePhiT = np.tan(i) * ( ( 1+ w*np.tan(THETA) )/( w-np.tan(THETA) ) )
 #       mask = np.abs(SenodePhiT) <= 1.  #Discard non real values for the angle of tangent line
        SenodePhiT[np.abs(SenodePhiT)>=1.] =np.nan #other way to set mask
        mask = (np.sin(PHI-0.5*DPHI) <= SenodePhiT)  &  (SenodePhiT <= np.sin(PHI +0.5*DPHI))
        xi,yi,zi = x*np.cos(i) - z*np.sin(i), y , x*np.sin(i)+z*np.cos(i) 
        plt.plot(xi[mask],yi[mask],'*')
        plt.axis([-5.0,0.5, -5.1, 5.1])
        plt.xlabel('z')
        plt.ylabel('r')
        plt.title('Bowshock shapes for {} wind and beta = {} (Test with i={})'.format(inn,b,i))
        plt.savefig('{}-wind-plane-sky-beta = {} (Test i= {} ).png'.format(inn,b,i))
        plt.clf()
