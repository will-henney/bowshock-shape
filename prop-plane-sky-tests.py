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
#parser = argparse.ArgumentParser(description="Rotation angle")
#parser.add_argument("--iangle", "-i", default=0.0, type=float, help="Rotation of z axis in radians")
#cmdargs = parser.parse_args()
#i = cmdargs.iangle

N = 200
#h=1e-6
theta = np.linspace(0,np.pi,N)
beta = np.array([0.15, 0.05, 0.02, 0.005,0.2,0.3])
innertype = ['isotropic','proplyd']
#phi = np.linspace(-np.pi,np.pi,N,endpoint = False) #Phi is the azimutal angle in spherical coordinates to make the revolution surface of proplyd
#Meshgrid theta and phi:
#THETA,PHI = np.meshgrid(theta,phi)
#DPHI = 2*np.pi/N
i = np.linspace(0.0,0.5236,8)
for inn in innertype:
    print "******{} case******".format(inn)
    for b in beta:
        for j in i:
            shell = Shell(beta=b,innertype=inn)
#           Calculing R(theta)
            R=shell.radius(theta)
#           Rdr = shell.radius(theta + h)
            R[R<=0] = np.nan # Set neagtive R to NaN
#           W = OMEGA(R,Rdr,h)
            w = omega(R,theta)
#           x,y,z = R*np.cos(THETA),R*np.sin(THETA)*np.cos(PHI),R*np.sin(THETA)*np.sin(PHI)
            SenodePhiT = np.tan(j) * ( ( 1+ w*np.tan(theta) )/( w-np.tan(theta) ) )
 #          mask = np.abs(SenodePhiT) <= 1.  #Discard non real values for the angle of tangent line
            SenodePhiT[np.abs(SenodePhiT)>=1.] =np.nan #other way to set mask
#           mask = (np.sin(PHI-0.5*DPHI) <= SenodePhiT)  &  (SenodePhiT <= np.sin(PHI +0.5*DPHI))
            xi,yi=R*np.cos(theta)*np.cos(j)-R*np.sin(theta)*SenodePhiT*np.sin(j),R*np.sin(theta)*np.sqrt(1-SenodePhiT**2) #Uses eq. A9 from Henney et al 
#                                                                           instead eq A1, A2. The mask in this case 
            plt.plot(xi,yi,label= 'i={}'.format(j))                             #becomes unnecessary
        plt.axis([-1.0,1.0,-0.05,2.0])
        plt.legend()
        plt.xlabel('z')
        plt.ylabel('r')
        plt.title('Bowshock shapes for {} wind and beta = {}'.format(inn,b))
        plt.savefig('{}-wind-plane-sky-beta-{}.png'.format(inn,b))
        plt.clf()
