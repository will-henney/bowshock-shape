from __future__ import print_function
from equation6 import Shell 
import numpy as np
import astropy
from astropy import coordinates as coord
from region_utils import region_point_to_string
import argparse
import json
from misc_utils import run_info
from scipy.optimize import bisect
import sys
import lmfit
"""
Goal: The idea of this program is create a CRW bowshock and write an output *xy.json file to use it later with Will's former programs
Steps: 
1.- Create the bowshock arc with the equation6 library for a fixed beta parameter for isotropic and/or proplyd case (perform to any beta value for both 
cases).

3.- Use the astropy packages as well other libraries to write the coordinates of the arc in a ds9 region format  
"""
def mirror(R,t):
    """
    Creating the other part of bowshock, assuming it is symetrical
    """
    
    tf = np.array([t,-t])
    rf = np.array([R,R])
    tf,rf = tf.reshape(2*len(t),), rf.reshape(2*len(R),)
 
    def find_th_order(t):
        """
        Find the correct order of theta
        """ 
        return t.argsort()
    order = find_th_order(tf)
    tf = tf[order]
    rf = rf[order]

    return rf,tf


def inner(r,t,m):
    """
    Computes the shell width and returns the inner shell shape
    """
    
    def omega(r,t):
        """
        Derivative of log(r) respect theta
        """
        ww = np.zeros_like(r)
        ww[:-1] = np.diff(np.log(r))/np.diff(t)
        ww[-1]=ww[-2]  # in the wings the slope is more or less constant (instead of zero)
        return ww

    def alpha(w,t):
        ta = (1+w*np.tan(t))/(w-np.tan(t))
        return -np.arctan(ta)

    H0 = 3./(4*m**2+7)
    n = 3.4*(4*m**2+1)/(4*m**2+35)
    psi = t+alpha(omega(r,t),t)-0.5*np.pi
    psi[m*np.cos(psi)<1.]=np.nan
    H = H0*np.cos(psi)**(-n)
    return r- r[0]*H/np.cos(psi)

def rotate(i,r,t):
    """
    shell shape rotated by an angle inc respect the "plane of sky"
    """
    x = r/np.cos(i)*np.cos(t)
    y = r/np.cos(i)*np.sin(t)
    

    def omega(r,t):
        """
        Derivative of log(r) respect theta
        """
        ww = np.zeros_like(r)
        ww[:-1] = np.diff(np.log(r))/np.diff(t)
        ww[-1]=ww[-2]  # in the wings the slope is more or less constant (instead of zero)
        return ww

    def sin_tan_phi(i,w,t):
        """
        sin of phi angle at the tangent line
        """
        return np.tan(i)*((1+w*np.tan(t))/(w-np.tan(t)))

    w = omega(r,t)
    s = sin_tan_phi(i,w,t)
    s[s>1.0] = np.nan
    xr,yr =  x*np.cos(i) - np.sin(t)*s*np.sin(i),y*np.sqrt(1-s**2) 
    mask = np.isfinite(xr) & np.isfinite(yr)
    return xr[mask],yr[mask]

def theta_lim(beta):
    "Asymptotic opening angle from CRW eq (28)"
    def f(theta):
        "Function to be zeroed: f(theta) = 0 for theta = theta_lim"
        return theta - np.tan(theta) - np.pi/(1.0 - beta)
    return bisect(f, 0.5*np.pi+0.01, np.pi)

def radius_from_point(x, y, x0, y0):
    return np.hypot(x-x0, y-y0)


def Rc_from_data(x, y, xc, yc):
    return np.mean(radius_from_point(x, y, xc, yc))


def deviation_from_circle(x, y, xc, yc):
    return radius_from_point(x, y, xc, yc) - Rc_from_data(x, y, xc, yc)


def model_minus_data(params, x, y):
    xc = params["xc"].value
    yc = params["yc"].value
    return deviation_from_circle(x, y, xc, yc)


def fit_circle(x, y, xc=0.0, yc=0.0):
    """
    Fit a circle to the shape of an arc with coordinates x, y

    Optionally provide initial guesses for the circle parameters: 
    xc, yc, Rc
    """
    params = lmfit.Parameters()
    params.add("xc", value=xc)
    params.add("yc", value=yc)
    lmfit.minimize(model_minus_data, params, args=(x, y))
    lmfit.report_errors(params)
    xc = params["xc"].value
    yc = params["yc"].value
    Rc = Rc_from_data(x, y, xc, yc)
    return Rc, xc, yc


parser = argparse.ArgumentParser(description = "model parameters")
parser.add_argument("--beta",type= float,default=0.01,help="two winds momentum ratio")
parser.add_argument("--type",type=str,default="isotropic",help="type of inner wind")
parser.add_argument("--mach",type=float,default=3.0,help="Mach number of inner wind")
cmd_args = parser.parse_args()

#1
beta = cmd_args.beta
if cmd_args.type=="isotropic":
    theta = np.linspace(0,theta_lim(beta),1500)
elif cmd_args.type=="proplyd":
    theta = np.linspace(0,np.pi,1500)
else:
    print("Unknown shell type")
    sys.exit()
shell = Shell(beta,cmd_args.type)
R = shell.radius(theta)
R[R<0]=np.nan
R_in = inner(R,theta,cmd_args.mach)
R_com,theta_com = mirror(R_in,theta)





arcdata = {
    "info":{
        "Description": "CRW theoretical stationary bowshocks",
        "history":["Initially created by "+ run_info()],
        "beta": beta,
        "inner wind type": cmd_args.type
        },
    "help":{
        "M": "Mach number of inner wind",
        "beta":"Two-winds momentum ratio",
        "i": "Rotation angle between proplyd reference frame and observer reference frame",
        "inner wind type": "Angular variation of density's inner wind",
        "x": "(list) x-coordinates of the outer shell normalized with distance",
        "y": "(list) y-coordinates of the outer shell normalized with distance"
        },
    "inner":{}
}
for inc in [0.0,15.0,30.0,45.0,60.0]:
    xp,yp = rotate(np.radians(inc),R_com,theta_com) 
    m = np.abs(np.degrees(np.arctan2(yp, xp))) - 45.0 <= 0.
    try:
        rc,xc,yc = fit_circle(xp[m],yp[m])
    except TypeError:
        break
    arcdata["inner"][str(inc)] = {
        "M":cmd_args.mach,
        "x":list(xp),
        "y":list(yp),
        "Rc": rc,
        "xc": xc,
        "yc": yc,
        "PAc":np.degrees(np.arctan2(-yc,-xc))        
    }      


jsonfile= "inner-{}-beta-{}-m-{}-arcdata.json".format(cmd_args.type,beta,cmd_args.mach)
with open(jsonfile,"w") as f:
    json.dump(arcdata,f,indent=4)


