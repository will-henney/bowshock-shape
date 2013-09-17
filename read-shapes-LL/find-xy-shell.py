import numpy as np
import argparse
import json

parser = argparse.ArgumentParser(
    description="""Find (X, Y) positions of shell boundaries from a DS9 region file""")

parser.add_argument("--region", type=str,
                    default="LL1-forma.reg",
                    help="Region file containing shell and star positions")

parser.add_argument("--debug", action="store_true",
                    help="Print out verbose debugging info about each line in region file")

cmd_args = parser.parse_args()
regionfile = cmd_args.region


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


# lists to contain x, y coords of inner and outer arcs
inner_x, inner_y = [], []
outer_x, outer_y = [], []
star_x, star_y = None, None
with open(regionfile) as f:
    lines = f.readlines()
    for line in lines:
        skipthisline = line.startswith("#") \
            or line.startswith("global") \
            or not "#" in line
        if skipthisline:
            continue
        ra, dec, point_type, text = extract_data(line)

        sdeg, samin, sasec = dec.split(':')
        deg, amin, asec = float(sdeg), float(samin), float(sasec)
        dec_arcsec = 3600.*deg + np.sign(deg)*(60*amin + asec)
        shr, smin, ssec = ra.split(':')
        hr, mn, sec = float(shr), float(smin), float(ssec)
        ra_arcsec = 15*np.cos(np.radians(dec_arcsec/3600.0))*(3600.*hr+60.*mn+sec)
        

        if point_type == "circle":
            # Position of star
            star_x, star_y = ra_arcsec, dec_arcsec
        elif point_type == "cross":
            # Points that trace inner edge
            inner_x.append(ra_arcsec)
            inner_y.append(dec_arcsec)
        elif point_type == "x":
            # Points that trace outer edge
            outer_x.append(ra_arcsec)
            outer_y.append(dec_arcsec)

assert star_x is not None, "Central star position not found - sorry!"

inner_x = np.array(inner_x) - star_x
inner_y = np.array(inner_y) - star_y
outer_x = np.array(outer_x) - star_x
outer_y = np.array(outer_y) - star_y

arc_data = {}

for arc_type, x, y in [
        ["inner", inner_x, inner_y],
        ["outer", outer_x, outer_y],
]:
    R = np.hypot(x, y)
    th = np.arctan2(y, x)
    # Need to make sure all arrays sorted in ascending theta order
    order = th.argsort()
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
    # Find closest point to star
    i0 = np.argmin(R)
    assert i0 > 1 and i0 + 1 < len(R), "Not enough points either side of closest point: i0 = {}".format(i0)
    # Fit quadratic through radii @ i0-1, i0, i0+1
    nbhood = slice(i0-1, i0+2)
    p = np.poly1d(np.polyfit(th[nbhood], R[nbhood], 2))
    # The theta that minimizes R is the (only) root of the derivative of p
    th0, = p.deriv().r
    R0 = p(th0)
    # Extract 2nd derivative
    d2R_dth2, = p.deriv().deriv()
    # And check that th0 is really a minimum of R(th)
    assert d2R_dth2 > 0.0, "Polynomial\n {:s}\n\ndoes not have a minimum!".format(p)
    # Transform to new frame where x-axis is along the th0 direction
    xx = R*np.cos(th - th0)
    yy = R*np.sin(th - th0)
    # Save results into data structure
    arc_data[arc_type] = {
        "theta0": np.degrees(th0),
        "R0": R0,
        "x": list(xx),
        "y": list(yy),
        }

jsonfile = regionfile.replace("-forma.reg", "-xy.json")
with open(jsonfile, "w") as f:
    json.dump(arc_data, f, indent=4)

