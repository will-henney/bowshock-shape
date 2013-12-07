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


def update_arc_data(data):
    """
    Update the data dictionary for an arc with the hyperbola fit parameters
    """
    x, y = np.array(data["x"]), np.array(data["y"])
    # Force the hyperbola axis to be along the circle fit axis
    # Initial guess is that center is also the same.
    Rh, thh, PAh, xh, yh = conic_utils.fit_hyperbola(
        x, y, Rh=data["Rc"], thh=45.0, PAh=data["PAc"], xxh=data["xc"], yyh=data["yc"])
    data.update(Rh=Rh, thh=thh, PAh=PAh, xh=xh, yh=yh)
    if cmd_args.savefig:
        plt.plot(-x, y, ".")
        print(arc, ":", xh, yh, Rh, thh)
        plt.plot(-xh, yh, "+" + colors[arc], ms=5.0)
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
infile = cmd_args.source + "-xyc.json"

db = json.load(open(infile))

for arc in "inner", "outer":
    update_arc_data(db[arc])

db["info"]["history"].append("Hyperbola fits added by " + run_info())

# save database to the same name: *-xyc.json
outfile = infile
with tempfile.NamedTemporaryFile(
        'w', dir=os.path.dirname(outfile), delete=False) as tf:
    json.dump(db, tf, indent=4)
    tempname = tf.name
os.rename(tempname, outfile)

if cmd_args.savefig:
    plotfile = infile.replace("-xy.json", "-conic-fits.pdf")
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

