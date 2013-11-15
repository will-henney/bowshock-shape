import numpy as np
import argparse
import json
from astropy import coordinates as coord
from astropy import units as u
from angle_utils import vector_separation
from misc_utils import run_info

parser = argparse.ArgumentParser(
    description="""Find (X, Y) positions of shell boundaries from a DS9 region file""")

parser.add_argument("source", type=str,
                    default="LL1",
                    help="Name of source, taken as prefix for region file containing shell and star positions")

parser.add_argument("--debug", action="store_true",
                    help="Print out verbose debugging info about each line in region file")

cmd_args = parser.parse_args()
regionfile = cmd_args.source + "-forma.reg"


def extract_data(line):
    """Find the coordinates, point type,
    and description text from the
    line of a ds9 regions file

    """

    coordstring, paramstring = line.split("#")
    shape, numbers = coordstring.split(")")[0].split("(")
    ra, dec = numbers.split(",")[0:2]
    shape = shape.strip()
    if cmd_args.debug:
        print "Shape: ({})".format(shape)
        print "RA, Dec: ", ra, dec
        print "Params: ", paramstring

    if "point" in shape:
        point_type = paramstring.split("point=")[1].split(" ")[0].strip()
    else:
        point_type = None

    if "tag" in line:
        text = paramstring.split("tag={")[1].split("}")[0].strip()
    elif "text" in line:
        text = paramstring.split("text={")[1].split("}")[0].strip()
    else:
        text = "NO TEXT"

    return ra, dec, point_type, text

th1C_RA_dg = 83.818289
th1C_Dec_dg = -5.3895909
th1C_c = coord.ICRSCoordinates(th1C_RA_dg, th1C_Dec_dg, unit=(u.degree, u.degree))

# lists to contain coords of inner and outer arcs
inner_c, outer_c = [], []
star_c = None



with open(regionfile) as f:
    lines = f.readlines()
    for line in lines:
        skipthisline = line.startswith("#") \
            or line.startswith("global") \
            or not "#" in line
        if skipthisline:
            continue
        ra, dec, point_type, text = extract_data(line)
        
        c = coord.ICRSCoordinates(ra, dec, unit=(u.hour, u.degree))

        if point_type == "circle":
            # Position of star
            star_c = c
            ra_star, dec_star = ra, dec
            D_star, pa_star = vector_separation(star_c, th1C_c, mode="polar")
            # PA must be in radians for next steps, then converted
            # back to degrees at the end
            pa_star = np.radians(pa_star)

        elif point_type == "cross":
            # Points that trace inner edge
            inner_c.append(c)
        elif point_type == "x":
            # Points that trace outer edge
            outer_c.append(c)

assert star_c is not None, "Central star position not found - sorry!"

inner_x, inner_y = [], []
for c in inner_c:
    x, y = vector_separation(star_c, c)
    inner_x.append(x)
    inner_y.append(y)
inner_x = np.array(inner_x)
inner_y = np.array(inner_y)

outer_x, outer_y = [], []
for c in outer_c:
    x, y = vector_separation(star_c, c)
    outer_x.append(x)
    outer_y.append(y)
outer_x = np.array(outer_x)
outer_y = np.array(outer_y)

arc_data = {
    "star": {
        "id": regionfile.replace("-forma.reg", ""),
        "RA": ra_star, 
        "Dec": dec_star,
        "RA_dg": star_c.ra.deg, 
        "Dec_dg": star_c.dec.deg,
        "PA": np.degrees(pa_star) % 360.0,
        "D": D_star
     } 
}


def find_th_order(th): 
    """Returns a sort order for a collection of angles theta
    
    Takes care to account for the wrap-around of angles by shifting
    the star-th1C vector to be at pi, so that all points are (with
    luck) in the range [0, 2 pi]

    """
    th1 = (canonicalize(th - pa_star) + np.pi) % (2*np.pi)
    if cmd_args.debug: 
        print "Finding theta order: " 
        print "    th = ", np.degrees(th)
        print "    pa_star = ", np.degrees(pa_star)
        print "    th1 = ", np.degrees(th1)
        print "    order = ", th1.argsort()
    return th1.argsort()


def canonicalize(th, unit="radians"):
    """Fold an angle theta into the canonical range [-pi:pi]"""
    if unit == "radians":
        return ((th + np.pi) % (2*np.pi)) - np.pi
    elif unit == "degrees":
        return ((th + 180.0) % (360.0)) - 180.0
    else:
        raise NotImplementedError

        
for arc_type, x, y in [
        ["inner", inner_x, inner_y],
        ["outer", outer_x, outer_y],
]:
    if len(x) == 0:
        continue
    R = np.hypot(x, y)
    th = np.arctan2(x, y) % (2*np.pi) # this now a PA for simplicity
    # Need to make sure all arrays sorted in ascending theta order
    
    order = find_th_order(th)
    x = x[order]
    y = y[order]
    th = th[order]
    R = R[order]
    if cmd_args.debug:
        print arc_type
        print "x: ", x
        print "y: ", y
        print "R: ", R
        print "th: ", np.degrees(th)

    i0 = np.argmin(R)
    if i0 > 0 and i0 + 1 < len(R):
        nbhood = slice(i0-1, i0+2)
    else:
        nbhood = slice(None)    # Use all the points in the quadratic fit
        print "Warning: Closest point of {} arc is at one end, using all points in parabola fit".format(arc_type)

    p = np.poly1d(np.polyfit(canonicalize(th[nbhood] - th[i0]), R[nbhood], 2))
    # The theta that minimizes R is the (only) root of the derivative of p
    th0, = p.deriv().r
    R0 = p(th0)
    # And check that th0 is really a minimum of R(th)
    assert p[0] > 0.0, "Polynomial\n {:s}\n\ndoes not have a minimum!".format(p)
    th0 += th[i0]

    # Transform to new frame where x-axis is along the th0 direction
    xx = R*np.cos(th - th0)
    yy = R*np.sin(th - th0)
    PA0 = np.degrees(th0) % 360.0

    if cmd_args.debug:
        print "R0 = {:.2f} arcsec, PA0 = {:.2f} deg".format(R0, PA0)
    # Save results into data structure
    arc_data[arc_type] = {
        "PA0": PA0,
        "R0": R0,
        "x": list(x),
        "y": list(y),
        "X": list(xx),
        "Y": list(yy),
        "R": list(R),
        "theta": list(np.degrees(canonicalize(th - th0))),
        }

arc_data["help"] = {
    "arcs":{
        "PA0": "[deg] Position angle of symmetry axis",
        "R0": "[arcsec] Radius along symmetry axis",
        "x": "(list) [arcsec] RA offsets from star of shell points",
        "y": "(list) [arcsec] Dec offsets from star of shell points",
        "X": "(list) [arcsec] Offsets from star parallel to axis of shell points",
        "Y": "(list) [arcsec] Offsets from star perpendicular to axis of shell points",
        "R": "(list) [arcsec] Radii of shell points from star",
        "theta": "(list) [deg] Angle from axis of shell points",
        "h": "[arcsec] Bowshock width along symmetry axis"
    },
    "star": {
        "id": "Source identifier deduced from region filename",
        "RA": "[HH:MM:SS] Right ascension of source", 
        "Dec": "[deg:arcmin:arcsec] Declination of source",
        "RA_dg": "[deg] Right ascension of source", 
        "Dec_dg": "[deg] Declination of source",
        "PA": "[deg] Position Angle of th^1 C from source",
        "D": "[arcsec] Distance of th^1 C from source",
     },
    "thickness": {
        "h0": "Shell thickness on axis: h0 = R0(outer) - R0(inner)"
    }

}

arc_data["info"] = {
    "description": "JSON data file for stationary bowshock arcs",
    "history": ["Initially created by " + run_info()],
    "region": regionfile,
} 

arc_data["thickness"] = {
    "h0": np.abs(arc_data["inner"]["R0"]-arc_data["outer"]["R0"]),
} 

jsonfile = regionfile.replace("-forma.reg", "-arcdata.json")
with open(jsonfile, "w") as f:
    json.dump(arc_data, f, indent=4)













