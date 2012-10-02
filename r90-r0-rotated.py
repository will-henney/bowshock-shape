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
i = np.linspace(-20,20,9)


N = 100
theta = np.linspace(0,0.5*np.pi,N)
beta = np.linspace(0.01,0.5,20)
innertype = ['isotropic','proplyd']
#phi = np.linspace(-np.pi,np.pi,N,endpoint = False) #Phi is the azimutal angle in spherical coordinates to make the revolution surface of proplyd
#Meshgrid theta and phi:
#THETA,PHI = np.meshgrid(theta,phi)
#DPHI = 2*np.pi/N
#i = np.linspace(-30.0,30.0,8)

for inn in innertype:
    print "******{} case******".format(inn)
    for j in i:
        R0 = list()
        R90 = list()
        for b in beta:
            shell = Shell(beta=b, innertype=inn)    
            R = shell.radius(theta)
            R[R<=0] = np.nan
            w = omega(R,theta)
#           x,y,z = R*np.cos(THETA),R*np.sin(THETA)*np.cos(PHI),R*np.sin(THETA)*np.sin(PHI)
            SenodePhiT = np.tan(j*np.pi/180.) * ( ( 1+ w*np.tan(theta) )/( w-np.tan(theta) ) )
            SenodePhiT[np.abs(SenodePhiT)>=1.] =np.nan #other way to set mask
#           mask = (np.sin(PHI-0.5*DPHI) <= SenodePhiT)  &  (SenodePhiT <= np.sin(PHI +0.5*DPHI))
#           xi,yi,zi = x*np.cos(j*np.pi/180.)-z*np.sin(j*np.pi/180.),y,x*np.sin(j*np.pi/180.)+z*np.cos(j*np.pi/180.)Convention sign in article Henney et al
            xi,yi=R*np.cos(theta)*np.cos(j*np.pi/180.)-R*np.sin(theta)*SenodePhiT*np.sin(j*np.pi/180.),R*np.sin(theta)*np.sqrt(1-SenodePhiT**2)
            xim,yim=xi[np.isnan(xi)==False],yi[np.isnan(yi)==False] #Removing nan elements from xi and yi
#           xim,yim = xi[mask],yi[mask] # x and y rotated and masked. 1D Arrays
            #print xim[0],yim[0]  Test #1
            #print xim[-1],yim[-1] Test#2
            R0.append(xim[0]) # Looks like is a secure criterium ( At least in the beta range (0,0.1] ). For high i and high beta (~0.5), odd things happen
            R90.append(yim[-1]) # This is even a more secure criterium
        #print "**********" 
        plt.plot(R0,R90,'o-',label = 'i={}'.format(j))
    plt.legend(loc="upper left")
    plt.axis([0,0.5,0,1.5])
    plt.xlabel("R_0 / D")
    plt.ylabel("R_90 / D")
    plt.title("Perpendicular versus parallel bowshock radii")
    plt.savefig("{}-shell-test-R0-R90.png".format(inn))
    plt.clf()
