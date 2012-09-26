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
#parser = argparse.ArgumentParser(description="Rotation angle")#
#parser.add_argument("--iangle", "-i", default=0.0, type=float, help="Rotation of z axis in radians")
#cmdargs = parser.parse_args()
#i = cmdargs.iangle
i = np.linspace(0,45,9)


N = 50
h=1e-6
theta = np.linspace(0,0.5*np.pi,N)
beta = np.linspace(0.,0.1,10)
innertype = ['isotropic','proplyd']
phi = np.linspace(-np.pi,np.pi,N,endpoint = False) #Phi is the azimutal angle in spherical coordinates to make the revolution surface of proplyd
#Meshgrid theta and phi:
THETA,PHI = np.meshgrid(theta,phi)
DPHI = 2*np.pi/N
i = np.linspace(0.0,30.0,5)

for inn in innertype:
    print "******{} case******".format(inn)
    for j in i:
        R0 = list()
        R90 = list()
        for b in beta:
            shell = Shell(beta=b, innertype=inn)    
            R = shell.radius(theta)
            w = omega(R,theta)
            x,y,z = R*np.cos(THETA),R*np.sin(THETA)*np.cos(PHI),R*np.sin(THETA)*np.sin(PHI)
            SenodePhiT = np.tan(j*np.pi/180.) * ( ( 1+ w*np.tan(THETA) )/( w-np.tan(THETA) ) )
            SenodePhiT[np.abs(SenodePhiT)>=1.] =np.nan #other way to set mask
            mask = (np.sin(PHI-0.5*DPHI) <= SenodePhiT)  &  (SenodePhiT <= np.sin(PHI +0.5*DPHI))
            xi,yi,zi = x*np.cos(j*np.pi/180.) - z*np.sin(j*np.pi/180.), y , x*np.sin(j*np.pi/180.)+z*np.cos(j*np.pi/180.) #Changed i sign
            xim,yim = xi[mask],yi[mask] # x and y rotated and masked. 1D Arrays
            R0.append(xim[0]) # R0 is the first element, it has only xi component 
            R90.append(yim[-1]) # R90 is the last element, it has only yi component
        plt.plot(R0,R90,'o-',label = 'i={}'.format(j))
    plt.legend()
    plt.axis([0,0.5,0,1.5])
    plt.xlabel("R_0 / D")
    plt.ylabel("R_90 / D")
    plt.title("Perpendicular versus parallel bowshock radii")
    plt.savefig("{}-shell-test-R0-R90.png".format(inn))
    plt.clf()
