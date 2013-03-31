import numpy as np
import matplotlib.pyplot as plt
from equation6 import Shell
import argparse
from scipy.optimize import fsolve,bisect,leastsq
from scipy.interpolate import interp1d

"""
In this script, we'll create the proyection of proplyd in the sky plane, assuming the plane of proplyd is rotated respect plane of sky by an angle i
"""

def calc_R(xc, yc):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((x_fit-xc)**2 + (y_fit-yc)**2)


def f_2(c):
    """ calculate the algebraic distance between the 2D points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(*c)
    return Ri - Ri.mean()

def theta_lim(beta): 
    "Asymptotic opening angle from CRW eq (28)"
    def f(theta):
	"Function to be zeroed: f(theta) = 0 for theta = theta_lim"
	return theta - np.tan(theta) - np.pi/(1.0 - beta)
    return bisect(f, 0.5*np.pi+0.01, np.pi)

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

parser = argparse.ArgumentParser(description = """Find the residual between radius of bowshock model and the circle fit as function of theta""")

parser.add_argument("--tdata",type=float,default = 45,help = 'upper angle of data for fitting')

cmd_args = parser.parse_args()

tdata = cmd_args.tdata

Nth = 200
Ninc = 10

beta = np.array([0.16, 0.08, 0.04, 0.02, 0.01])
innertype = ['isotropic','proplyd']
lw = dict(isotropic = 1, proplyd = 2)
color = {'isotropic':'b.','proplyd':'r.'}
incdeg = np.array([0.0, 15.0, 30.0, 45.0, 60.0, 75.0])
colors = "bgrmkcy"

for b in beta:
    print "******beta = {}******".format(b)
    thlim = theta_lim(b) 
    theta = np.linspace(0,thlim,Nth)
    for inn in innertype:
        inc = np.radians(incdeg)
        shell = Shell(beta=b,innertype=inn)
#       Calculing R(theta)
        R=shell.radius(theta)
        R[R<=0] = np.nan # Set neagtive R to NaN
        w = omega(R,theta)
#        for jnc, col in zip(inc, colors[:len(inc)]):
        SenodePhiT = np.tan(inc[0]) * ( ( 1+ w*np.tan(theta) )/( w-np.tan(theta) ) )
        SenodePhiT[np.abs(SenodePhiT)>=1.] =np.nan 
        xi = (R/np.cos(inc[0]))*(np.cos(theta)*np.cos(inc[0]) 
                              - np.sin(theta)*SenodePhiT*np.sin(inc[0]))
        yi = (R/np.cos(inc[0]))*np.sin(theta)*np.sqrt(1-SenodePhiT**2) 
        mask = np.isfinite(xi) & np.isfinite(yi)
        xim, yim = xi[mask], yi[mask]
        xplot,yplot = np.array([xim,xim]),np.array([-yim,yim])
	x2,y2 = xplot.reshape(2*len(xim),),yplot.reshape(2*len(xim),)
        if inn == "isotropic":
            label = 'i={:.0f}'.format(np.degrees(inc[0]))
        else:
            label = None
        x_m = 0.
        y_m = 0.
        center_estimate = x_m, y_m
        x_fit, y_fit=x2[(np.abs(np.degrees(np.arctan2(y2,x2)))) -tdata <= 0.], y2[(np.abs(np.degrees(np.arctan2(y2,x2)))) -tdata <= 0.]
        center_2, ier = leastsq(f_2, center_estimate)
        xc_2, yc_2 = center_2
        Ri_2       = calc_R(xc_2, yc_2)
        R_2        = Ri_2.mean()
       
        theta_fit = np.linspace(-np.pi, np.pi,500)
        x_fit2 = xc_2 + R_2*np.cos(theta_fit)
        y_fit2 = yc_2 + R_2*np.sin(theta_fit)
	#calculating the residuals

	#Shock data
	t_m = np.arctan2(y2[x2>=0.],x2[x2>=0.])
	RM = np.sqrt(x2[x2>=0.]**2 + y2[x2>=0.]**2)

	#circle data
	RC = np.sqrt(x_fit2[x_fit2>= (xc_2 + R_2*np.cos(np.radians(100.)))]**2 + y_fit2[x_fit2>=(xc_2 + R_2*np.cos(np.radians(100.)))]**2)
	t_c = np.arctan2(y_fit2[x_fit2 >=(xc_2 + R_2*np.cos(np.radians(100.)))],x_fit2[x_fit2>=(xc_2 + R_2*np.cos(np.radians(100.)))])
#        print np.degrees(t_m.min()),np.degrees(t_m.max()),np.degrees(t_c.min()),np.degrees(t_c.max())
	#interpolating stuff
        f = interp1d(t_c,RC)
	#array of RCs with shock's angles
        RC2 = f(t_m)
	res = (RM - RC2)/RM  
	plt.plot(np.degrees(t_m),res,color[inn],label = 'residuo-{}'.format(inn))
    plt.xlabel('theta')
    plt.ylabel('e')
    plt.legend(loc ='best')
    plt.grid()
    plt.title('Residual R vs theta-b-{}-w-{}'.format(b,tdata))
    plt.savefig('Residual-b-{}-w-{}.png'.format(b,tdata))
    plt.clf()

