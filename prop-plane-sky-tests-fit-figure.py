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
#Choose between normalizing with R0 or normalizing with D
parser = argparse.ArgumentParser(description='Choose normalized bowshock with R0')
parser.add_argument('--axis',type=str, choices = ('r/r0','r/D'),default='r/D',help='Normalizing bowshock')
cmd_args = parser.parse_args()
isNormalized = cmd_args.axis == 'r/r0'

Nth = 200
Ninc = 10

beta = np.array([0.16, 0.08, 0.04, 0.02, 0.01])
innertype = ['isotropic','proplyd']
lw = dict(isotropic = 1, proplyd = 2)

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
#Adding the lower part of bowshock
        xplot,yplot = np.array([xim,xim]),np.array([-yim,yim])
	x2,y2 = xplot.reshape(2*len(xim),),yplot.reshape(2*len(xim),)
        xr0,yr0 = np.array([0.0,x2.max()]),np.array([0.0,y2[x2==x2.max()][-1]])
        X = x2/xr0[-1] if isNormalized else x2
        Y = y2/xr0[-1] if isNormalized else y2
#Label for the shock curve
        if inn == "isotropic":
            label = 'i={:.0f}-{} shock'.format(np.degrees(inc[0]),inn)
        else:
            label = 'i={:.0f}-{} shock'.format(np.degrees(inc[0]),inn)
#Fitting the curve model
        x_m = 0.
        y_m = 0.
        center_estimate = x_m, y_m
        x_fit, y_fit=x2[(np.abs(np.degrees(np.arctan2(y2,x2)))) -45. <= 0.], y2[(np.abs(np.degrees(np.arctan2(y2,x2)))) -45. <= 0.]
        center_2, ier = leastsq(f_2, center_estimate)
        xc_2, yc_2 = center_2
        Ri_2       = calc_R(xc_2, yc_2)
        R_2        = Ri_2.mean()
       
        theta_fit = np.linspace(-np.pi, np.pi, 180)
        x_fit2 = xc_2 + R_2*np.cos(theta_fit)
        y_fit2 = yc_2 + R_2*np.sin(theta_fit)
#Normalizing the curve with R0
        X_f = x_fit2/xr0[-1] if isNormalized else x_fit2
        Y_f = y_fit2/xr0[-1] if isNormalized else y_fit2
#Plotting the bowshock, the curvature radius and adding fancy labels in each case   
        if inn == 'isotropic':
            xrc,yrc = np.array([xc_2,xc_2 - 0.5*R_2]), np.array([yc_2,yc_2 + np.sqrt(3.)*R_2/2.])
            xrcn = xrc/xr0[-1] if isNormalized else xrc
            yrcn = yrc/xr0[-1] if isNormalized else yrc
	    plt.plot(X, Y,'b.', linewidth=lw[inn], label=label)  
            plt.plot(X_f, Y_f, 'k--', label='fit for isotropic shock', lw=2)
	    plt.plot(xrcn,yrcn,'y-',label='Rc-isotropic',linewidth = 2)
            plt.annotate("Rc-isotropic",xy=(0.5*(xrcn[0] + xrcn[-1]),0.5*(yrcn[0] + yrcn[-1])), xytext=(-20,20),fontsize='x-small',
                alpha=1.0, textcoords='offset points',ha='left',va='top',
                bbox= dict(boxstyle='round,pad=0.5',fc='green',alpha=0.5),
                arrowprops= dict(arrowstyle='->',connectionstyle='arc3,rad=0'))
        else:
            xrc,yrc = np.array([xc_2,xc_2 + R_2/np.sqrt(2.)]), np.array([yc_2,yc_2 + R_2/np.sqrt(2.)])
            xrcn = xrc/xr0[-1] if isNormalized else xrc
            yrcn = yrc/xr0[-1] if isNormalized else yrc
	    plt.plot(X, Y,'g.', linewidth=lw[inn], label=label)  
            plt.plot(X_f, Y_f, 'r--', label='fit for proplyd shock', lw=2)
            plt.plot(xrcn,yrcn,'c-',label = 'Rc-proplyd',linewidth=2)
            plt.annotate("Rc-proplyd",xy=(0.5*(xrcn[0] + xrcn[-1]),0.5*(yrcn[0] + yrcn[-1])),xytext=(20,20),fontsize='x-small',
                alpha=1.0, textcoords='offset points',ha='right',va='top',
                bbox= dict(boxstyle='round,pad=0.5',fc='green',alpha=0.5),
                arrowprops= dict(arrowstyle='->',connectionstyle='arc3,rad=0'))
#plotting the propld positions, the center of the circle and R0, with its fancy label
        plt.plot([xc_2/xr0[-1]], [yc_2/xr0[-1]], 'xg') if isNormalized else plt.plot([xc_2], [yc_2], 'xg')
        plt.plot(0.0,0.0,'xr')
	plt.plot([0.,1.],[0.,0.],'m-',linewidth = 2.) if isNormalized else plt.plot(xr0,yr0,'m-',linewidth = 2.)
        plt.annotate("R0",xy=(0.5*xr0[-1],0.5*yr0[-1]),xytext=(20,-20),fontsize='x-small',
            alpha=1.0, textcoords='offset points',ha='right',va='bottom',
            bbox= dict(boxstyle='round,pad=0.5',fc='green',alpha=0.5),
            arrowprops= dict(arrowstyle='->',connectionstyle='arc3,rad=0'))
#    plt.axis([-1.0,1.0,-0.05,yi[-1]+0.5])
    plt.legend(loc='best',prop=dict(size="x-small"))
    plt.grid()
    plt.axis("equal")
    if isNormalized:
        plt.xlim(-5.0, 2.0)
        plt.ylim(-5., 5.)
        plt.xlabel("z'/R0'")
        plt.ylabel("r'/R0'")
        plt.title('Normalized Bowshock shapes for wind and beta = {} and circle fit'.format(b))
        plt.savefig('Bowshock-wind-plane-sky-beta-{}-fitNorm.pdf'.format(b))
    else:
        plt.xlim(-1.5, 1.0)
        plt.ylim(-1., 1.)
        plt.xlabel("z'/D'")
        plt.ylabel("r'/D'")
        plt.title('Bowshock shapes for wind and beta = {} and circle fit'.format(b))
        plt.savefig('Bowshock-wind-plane-sky-beta-{}-fit3.pdf'.format(b))
    plt.clf()
