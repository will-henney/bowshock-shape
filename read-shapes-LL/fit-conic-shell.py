"""
Fit hyperbolae to bowshock arcs
"""
from __future__ import print_function

import numpy as np
import json
import argparse
import tempfile
import os
import matplotlib.pyplot as plt
from misc_utils import run_info
import conic_utils 


def update_arc_data(data,teo=False):
    """
    Update the data dictionary for an arc with the hyperbola fit parameters
    """
    x, y = np.array(data["x"]), np.array(data["y"])
    m = np.abs(np.degrees(np.arctan2(y, x))) - 100.0 <= 0.
    # Force the hyperbola axis to be along the circle fit axis
    # Initial guess is that center is also the same.
    if teo:
        Rh, thh, PAh, xh, yh = conic_utils.fit_conic(
            x[m], y[m], Rh=data["Rc"], thh=-45.0, PAh=data["PAc"], xxh=data["xc"], yyh=data["yc"])
    else:
        Rh, thh, PAh, xh, yh = conic_utils.fit_conic(
            x, y, Rh=data["Rc"], thh=45.0, PAh=data["PAc"], xxh=data["xc"], yyh=data["yc"])
    
    R0 = np.max(x)
    R90 = Rh*conic_utils.x90(R0,thh)

    data.update(Rh=Rh, thh=thh, PAh=PAh, xh=xh, yh=yh,R0=R0,R90=R90)
    if cmd_args.savefig:
        plt.plot(-x, y, ".")
        xx,yy = conic_utils.world_hyperbola(Rh,thh,PAh,xh,yh) 
        if teo:
            ax.grid()
            ax.set_xlim(-1.0,2)
            ax.set_ylim(-2,2)
            t=np.linspace(0,2*np.pi)
            ax.plot(-data["xc"], data["yc"], "+r", ms=5.0)
            ax.plot(-data["xc"]+data["Rc"]*np.cos(t),data["yc"]+data["Rc"]*np.sin(t),"g-",alpha=0.9)
            ax.plot(-xx,yy,"m-",alpha=1.0)
            ax.plot(xh,yh,"k+",ms=5.0)
        else:
            plt.plot(-xh, yh, "+" + colors[arc], ms=5.0)
            plt.plot(-xx,yy,colors[arc],alpha=0.9)
            print(arc, ":", xh, yh, Rh, thh)
        # c = plt.Circle((-xc, yc), radius=Rc,
        #                fc='none', ec=colors[arc], alpha=0.4, lw=0.2)
        # plt.gca().add_patch(c)
    return None


def PA_circle(xc, yc):
    """
    PA of star with respect to center of circle
    """
    return np.degrees(np.arctan2(-xc, -yc))


colors = {"inner": "m", "outer": "g"}

parser = argparse.ArgumentParser(
    description="""Fit conic sections to all the arcs and save to database file""")

parser.add_argument("--source", type=str,
                    default="LL1",
                    help="Name of source (prefix for data file)")
parser.add_argument("--savefig", action="store_true",
                    help="Save a figure showing the fit")
parser.add_argument("--debug", action="store_true",
                    help="Print out verbose debugging info")

cmd_args = parser.parse_args()
infile = cmd_args.source + "-arcdata.json"

db = json.load(open(infile))

if "beta" in infile:
    f = plt.figure()
    for i,inc in enumerate([0.0,15.0,30.0,45.0,60.0]):
        try:
            ax = f.add_subplot(2,3,i+1)
            update_arc_data(db["outer"]["i="+str(inc)],teo=True)
        except KeyError:
            #print "inclination not available",inc
            break
else:
    for arc in "inner", "outer":
        update_arc_data(db[arc])

db["info"]["history"].append("Conic fits added by " + run_info())

# save database to the same name: *-arcdata.json
outfile = infile
with tempfile.NamedTemporaryFile(
        'w', dir=os.path.dirname(outfile), delete=False) as tf:
    json.dump(db, tf, indent=4)
    tempname = tf.name
os.rename(tempname, outfile)

if cmd_args.savefig:
    plotfile = infile.replace("-arcdata.json", "-conic-fits.pdf")
    if "beta" in infile:
        f.tight_layout()
        f.set_size_inches(2*5,3*5)
        f.savefig(plotfile)
    else:
        plt.plot(0.0, 0.0, 'o')
        plt.axis("equal")
        plt.savefig(plotfile)



#
# Template for more robust file update
# Copied from http://blog.gocept.com/2013/07/15/reliable-file-updates-with-python/
#
# with tempfile.NamedTemporaryFile(
#       'w', dir=os.path.dirname(filename), delete=False) as tf:
#    tf.write(model.output())
#    tempname = tf.name
# os.rename(tempname, filename)

