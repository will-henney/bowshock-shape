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

parser = argparse.ArgumentParser(description = "model parameters")
parser.add_argument("--beta",type= float,default=0.01,help="two winds momentum ratio")
parser.add_argument("--type",type=str,default="isotropic",help="type of inner wind")
cmd_args = parser.parse_args()

#1
beta = cmd_args.beta
if cmd_args.type=="isotropic":
    theta = np.linspace(0,theta_lim(beta),1500)
elif cmd_args.type=="proplyd":
    theta = np.linspace(0,np.pi,1500)
else:
    print "Unknown shell type"
    sys.exit()
shell = Shell(beta,cmd_args.type)
R = shell.radius(theta)

R_com,theta_com = mirror(R,theta)
R[R<0]=np.nan



arcdata = {
    "info":{
        "Description": "CRW theoretical stationary bowshocks",
        "history":["Initially created by "+ run_info()],
        "beta": beta,
        "inner wind type": cmd_args.type
        },
    "help":{
        "beta":"Two-winds momentum ratio",
        "i": "Rotation angle between proplyd reference frame and observer reference frame",
        "inner wind type": "Angular variation of density's inner wind",
        "x": "(list) x-coordinates of the outer shell normalized with distance",
        "y": "(list) y-coordinates of the outer shell normalized with distance"
        },
    "outer":{}
}
for inc in [0.0,15.0,30.0,45.0,60.0,75.0]:
    xp,yp = rotate(np.radians(inc),R_com,theta_com) 
    arcdata["outer"]["i="+str(inc)] = {
        "x":list(xp),
        "y":list(yp)        
    }      


jsonfile= "{}-beta-{}-arcdata.json".format(cmd_args.type,beta)
with open(jsonfile,"w") as f:
    json.dump(arcdata,f,indent=4)


