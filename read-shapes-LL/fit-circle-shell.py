"""
Fit circles to bowshock arcs
"""


import numpy as np
import json
import argparse
import lmfit
from region_utils import region_hdr_lines, region_circle_to_string, region_point_to_string
from astropy import coordinates as coord
from astropy import units as u
import matplotlib.pyplot as plt

def create_arc_regions(arcdata):
    regions = []
    ra0 = coord.Longitude(arcdata["star"]["RA"], unit=u.hour)
    dec0 = coord.Latitude(arcdata["star"]["Dec"], unit=u.deg)
    for arc, c in ["inner", "magenta"], ["outer", "green"]:
        ra = ra0 + arcdata[arc]["xc"]*u.arcsec/np.cos(dec0.to(u.rad).value)
        dec = dec0 + arcdata[arc]["yc"]*u.arcsec
        radius = arcdata[arc]["Rc"]
        regions.append(region_circle_to_string(ra.to_string(sep=":"), dec.to_string(sep=":"), radius, text="", color=c))
        regions.append(region_point_to_string(ra.to_string(sep=":"), dec.to_string(sep=":"), "diamond", text=arc, color=c))
    return region_hdr_lines + regions


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


def update_arc_data(data):
    """
    Update the data dictionary for an arc with the circle fit parameters
    """
    x, y, theta = np.array(data["x"]), np.array(data["y"]), np.array(data["theta"])
    m = np.abs(theta) < 90.0
    if m.sum() < 3:
        print "Warning: only {} points within 90 deg of axis. Using all points instead".format(m.sum())
        m = np.ones_like(theta, dtype=bool)
    xc0 = data["R0"]*np.sin(np.radians(data["PA0"]+180))
    yc0 = data["R0"]*np.cos(np.radians(data["PA0"]+180))
    Rc, xc, yc = fit_circle(x[m], y[m], xc=xc0, yc=yc0)
    data.update(Rc=Rc, xc=xc, yc=yc, PAc=PA_circle(xc, yc))
    if cmd_args.savefig:
        plt.plot(-x, y, ".")
        plt.plot(-xc, yc, "+" + colors[arc], ms=5.0)
        c = plt.Circle((-xc, yc), radius=Rc,
                       fc='none', ec=colors[arc], alpha=0.4, lw=0.2)
        plt.gca().add_patch(c)
    return None


def PA_circle(xc, yc):
    """
    PA of star with respect to center of circle
    """
    return np.degrees(np.arctan2(-xc, -yc))


colors = {"inner": "m", "outer": "g"}

parser = argparse.ArgumentParser(
    description="""Fit circles to all the arcs and save as ds9 region file""")

parser.add_argument("infile", type=str,
                    default="LL1-xy.json",
                    help="JSON file containing arc data")
parser.add_argument("--savefig", action="store_true",
                    help="Save a figure showing the fit")
parser.add_argument("--debug", action="store_true",
                    help="Print out verbose debugging info")

cmd_args = parser.parse_args()
infile = cmd_args.infile

db = json.load(open(infile))

for arc in "inner", "outer":
    update_arc_data(db[arc])

outfile = infile.replace("-xy", "-xyc")
with open(outfile, "w") as f:
    json.dump(db, f, indent=4)

region_file = infile.replace("-xy.json", "-arcfits.reg")
with open(region_file, "w") as f:
      f.writelines([s + "\n" for s in create_arc_regions(db)])

if cmd_args.savefig:
    plotfile = infile.replace("-xy.json", "-arcfits.pdf")
    plt.plot(0.0, 0.0, 'o')
    plt.axis("equal")
    plt.savefig(plotfile)



