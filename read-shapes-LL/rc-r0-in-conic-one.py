"""
In this script, we create the projection of proplyd in the sky
plane, assuming the proplyd axis is at an angle i with respect to
plane of sky
"""

import argparse
import json
import conic_utils
import numpy as np
from scipy.optimize import bisect
from scipy.interpolate import interp1d
import lmfit
from equation6 import Shell
import matplotlib.pyplot as plt

#
# Command line arguments
#
parser = argparse.ArgumentParser(
    description="""Find parallel and conic fit for projected 
    bowshock models""")
parser.add_argument(
    '--tfit', type=float, default=90,
    help='upper angle of fit data')
parser.add_argument("--mach","-m",type=float,default=3.0,help="mach number of inner wind")
parser.add_argument("--beta",type=float,default=0.01,help = "winds momentum ratio")
parser.add_argument("--inc",type=float,default=0.0,help="inclination of plane of sky frame")
parser.add_argument("--innertype",type=str,default="isotropic",help="inner wind model")
cmd_args = parser.parse_args()


def theta_lim(beta):
    "Asymptotic opening angle from CRW eq (28)"
    def f(theta):
        "Function to be zeroed: f(theta) = 0 for theta = theta_lim"
        return theta - np.tan(theta) - np.pi/(1.0 - beta)
    return bisect(f, 0.5*np.pi+0.01, np.pi)


def omega(r, t):
    """
    Will's version of omega(theta), looks like work better.  uses
    d(log(r))/dtheta = 1/r * dr/dtheta (ingenious)
    """
    womega = np.zeros_like(r)
    womega[:-1] = np.diff(np.log(r))/np.diff(theta)
    return womega





#First, design a set of beta values to calculate R(theta) (Well, first
#than that, desgign a theta array)


beta = cmd_args.beta
Nth = 800
Ninc = 100
innertype = cmd_args.innertype
inc = cmd_args.inc

tfit = cmd_args.tfit
shelldata = {}
#Given the mach number, we can calculate the index of \cos\psi in the shell thickness and
#the thickness at the axis
H0 = 3./(4*cmd_args.mach**2+1.)
n  = 3.4*(4.*cmd_args.mach**2+1)/(4.*cmd_args.mach**2+35) 
       
thlim = theta_lim(beta)  
                             
theta = np.linspace(0, thlim, Nth)


# Initialize the data arrays to NaNs. R0 and R90 are enough to calculate
# the a,b parameters of an ellipse or hyperbola
R0 = np.nan
R90 = np.nan
thinf = np.nan
RE = np.nan
THINF = np.nan
x90 = np.nan

shell = Shell(beta=beta, innertype=innertype)
R = shell.radius(theta)
R[R<R[0]],theta[R<R[0]]=np.nan,np.nan
w = omega(R, theta)
tan_alpha = (1 + w*np.tan(theta)) /   (w - np.tan(theta))
alpha = -np.arctan(tan_alpha)
psi = theta+alpha-0.5*np.pi
psi_l = np.arccos(1./cmd_args.mach)
psi[psi>psi_l]=np.nan # For psi>psi_l The normal component of the inner wind velocity 
H = H0*np.cos(psi)**(-n) # is below the sound speed, so for those angles there is no shock
Rin = R- H*R[0]/np.cos(psi)
w_in = omega(Rin,theta)
tan_alpha_in = (1+w_in*np.tan(theta))/(w_in-np.tan(theta))
SenodePhiT = np.tan(inc)*tan_alpha_in 
# other way to set mask
SenodePhiT[np.abs(SenodePhiT) >= 1.] = np.nan

# Correct for projection and for fact that observed radii are normalised by
# D' = D cos(inc) (Now for inner radius)
xi = (Rin/np.cos(inc))*(np.cos(theta)*np.cos(inc)
     - np.sin(theta)*SenodePhiT*np.sin(inc))
yi = (Rin/np.cos(inc))*np.sin(theta)*np.sqrt(1-SenodePhiT**2)
mask = np.isfinite(xi) & np.isfinite(yi)
xim, yim = xi[mask], yi[mask]  # Removing nan elements from xi
                               # and yi

#returning to R' and theta' for calculating R90
Rp = np.hypot(xim,yim)
thetap = np.arctan2(yim,xim)
#Characterictic radii of shell
try:
    R0 = xim[0]
except IndexError:
    print("max inclination={}".format(inc))
try:
    f=interp1d(thetap,Rp,bounds_error=False)
    R90= f(0.5*np.pi)
except ValueError:
    print("max inclination={}".format(inc))

#Ellipse/hyperbola parameters got from R0 and R90 
EH = R90 < 2*R0 



if EH: #Ellipse case
    a = R0**2/(2*R0-R90)
    b = np.sqrt(a*R90)
    t=np.linspace(-0.5*np.pi,0.5*np.pi)
    xx,yy = a*np.cos(t) -np.sqrt(a**2-b**2), b*np.sin(t)
    thinf = np.degrees(np.arctan2(b,a))
elif 2*R0 == R90: #Parabola case
    xx,yy = R0*(1-2*t/np.pi),np.sqrt(4*R90*R0*t/np.pi)
    thinf = 0.0
else: # Hyperbola case
    a = np.abs(R0**2/(R90-2*R0)) 
    b = np.sqrt(a*R90)
    t=np.linspace(-np.pi,np.pi)
    xx,yy = -a*np.cosh(t) + np.hypot(a,b),-b*np.sinh(t)
    thinf = -np.degrees(np.arctan2(b,a))
# measure how similar are the conic and the bowshock
rr,tt = np.hypot(xx,yy),np.arctan2(yy,xx)
f1 = interp1d(thetap,Rp,bounds_error=False)
Rsh = f1(tt)

epsilon = (Rsh-rr)/Rsh

#Complete the bowshock to do an adequate fit with conic_utils
x2,y2 = np.array([xim,xim]),np.array([-yim,yim]) 
x2,y2 = x2.reshape(2*len(xim),),y2.reshape(2*len(yim),)
order = y2.argsort()
x2,y2 = x2[order],y2[order]
#Do the fitting with conic_utils
mask = np.abs(np.degrees(np.arctan2(y2,x2)))-cmd_args.tfit <= 0.0

#RE, THINF, PA,xc,yc = conic_utils.fit_conic(x2[mask],y2[mask],Rh=R0,thh=0.0,PAh=90.0,xxh=0.0,yyh=0.0)
#x90 = conic_utils.x90(R0,THINF)
#xx,yy = conic_utils.world_hyperbola(RE,THINF,PA,xc,yc)       
#plot outer shell, inner shell, the fit and characteristic radii

plt.rc("text", usetex=True)
plt.rc("font", family="serif")
f = plt.figure()
ax = f.add_subplot(2,1,1)
ax1 = f.add_subplot(2,1,2)
ax.plot(x2[mask],y2[mask],"r-",lw=2,label="M={}".format(cmd_args.mach))

#rotate and plot the outer shell
tan_alpha = (1+w*np.tan(theta))/(w-np.tan(theta)) 
SPT = np.tan(inc)*tan_alpha
x = (R/np.cos(inc))*(np.cos(theta)*np.cos(inc)
     - np.sin(theta)*SPT*np.sin(inc))
y = (R/np.cos(inc))*np.sin(theta)*np.sqrt(1-SPT**2)
x3,y3 = np.array([x,x]),np.array([-y,y]) 
x3,y3 = x3.reshape(2.*len(x),),y3.reshape(2.*len(y),)
order = y3.argsort()
x3,y3 = x3[order],y3[order]
ax.plot(x3,y3,"k-",lw=2.5,label=r"$\beta={}$ {}".format(beta,innertype))
#plt.plot(xx,yy,"c-",label = r"fit $\theta_{\infty}=$"+str(np.degrees(THINF)))
ax.plot([0.0,R0],[0.0,0.0],"b-",label=r"$R_0={:.3f}$".format(R0))
ax.plot(xx,yy,"m-",alpha=0.7,label=r"Conic with same $R_0$ and $R_{90}$ as shell")
ax.plot([0.0,0.0],[0.0,R90],"g-",label=r"$R_{90}=$"+str(R90))
#plt.plot([0.0,0.0],[0.0,x90],"y.-",alpha=0.7,label = r"fit $R_{90}=$"+str(x90))

ax.legend(loc="best",prop=dict(size="x-small"))
ax.grid()
ax.set_xlabel("x'/D'")
ax.set_ylabel("y'/D'")
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
ax.set_title(r"Shell shape. $\theta_{\infty}=$"+str(thinf))

ax1.plot(np.degrees(tt),epsilon,"r-")
ax1.grid()
ax1.set_xlabel(r"$\theta'$ (deg)")
ax1.set_ylabel(r"$\epsilon$")
ax1.set_xlim(0,np.degrees(thlim))
ax1.set_ylim(-0.05,0.05)
ax1.set_title("Difference between shock and conic")
f.set_size_inches(5,10)
f.savefig("shell-M={}-b={}-inc{}-type-{}.pdf".format(cmd_args.mach,beta,inc,innertype))
