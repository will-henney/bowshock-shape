"""
Fit circles to bowshock arcs
"""


import numpy as np
import json
import sys
import lmfit

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
    Rc, xc, yc = fit_circle(x[m], y[m], xc=x[m].mean(), yc=y[m].mean())
    data.update(Rc=Rc, xc=xc, yc=yc, PAc=PA_circle(xc, yc))
    return None

def PA_circle(xc, yc):
    """
    PA of star with respect to center of circle
    """
    return np.degrees(np.arctan2(-xc, -yc))


try: 
    infile = sys.argv[1]
except: 
    print "Usage: ", sys.argv[0], " INPUTFILE.json"
    sys.exit()

with open(infile) as f:
    db = json.load(f)


for arc in "inner", "outer":
    update_arc_data(db[arc])

outfile = infile.replace("-xy", "-xyc")

with open(outfile, "w") as f:
    json.dump(db, f, indent=4)






