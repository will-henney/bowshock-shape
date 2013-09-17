import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy import optimize

parser = argparse.ArgumentParser(
    description="""Find (X, Y) positions of shell boundaries from a DS9 region file""")

parser.add_argument("--region", type=str,
                    default="LL1-forma.reg",
                    help="Region file containing shell and star positions")

cmd_args = parser.parse_args()
regionfile = cmd_args.region
name, ext = regionfile.split('.')

def calc_R(xc, yc):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((x-xc)**2 + (y-yc)**2)


#@countcalls
def f_2(c):
    """calculate the algebraic distance between the 2D points and the
    mean circle centered at c=(xc, yc)"""
    Ri = calc_R(*c)
    return Ri - Ri.mean()


def extract_data(line):
    """Find the coordinates, point type,
    and description text from the
    line of a ds9 regions file

    """
    coordstring, paramstring = line.split("#")
    shape, numbers = coordstring.split(")")[0].split("(")
    ra, dec = numbers.split(",")[0:2]
    if shape == "point":
        point_type = paramstring.split("point=")[1].split(" ")[0]
    if "tag" in line:
        text = paramstring.split("tag={")[1].split("}")[0]
    elif "text" in line:
        text = paramstring.split("text={")[1].split("}")[0]
    else:
        text = "NO TEXT"
    return ra, dec, point_type, text


with open(regionfile) as f:
    lines = f.readlines()
    for line in lines:
        skipthisline = line.startswith("#") \
            or line.startswith("global") \
            or not "#" in line
        if skipthisline:
            continue
        ra, dec, point_type, text = extract_data(line)

        shr, smin, ssec = ra.split(':')
        hr, mn, sec = float(shr), float(smin), float(ssec)
        ra_arcsec = 15*(3600.*hr+60.*mn+sec)
        sdeg, samin, sasec = dec.split(':')
        deg, amin, asec = float(sdeg), float(samin), float(sasec)
        dec_arcsec = 3600.*deg + np.sign(deg)*(60*amin + asec)

        if point_type == "circle":
            # Position of star
        elif point_type == "cross":
            # Inner edge
        elif point_type == "x":
            # Outer edge

        if "shock" in text:
            if source not in Shapes:
                Shapes[source] = []
            Shapes[source].append(np.array([ra_arcsec, dec_arcsec]))
        else:
            Centers[source] = np.array([ra_arcsec, dec_arcsec])


# print Centers
# print
# print Shapes

proplyds = Shapes.keys()
proplyd = cmd_args.proplyd


# vector star->proplyd
try:
    vecD = Centers['th1C'] - Centers[proplyd]
except KeyError:
    # no center for this proplyd, so skip it
    print 'No center'
# scalar distance star->proplyd
D = np.hypot(*vecD)
# unit vector in direction star->proplyd
xunit = vecD/D
# unit vector perpendicular to star->proplyd
yunit = np.array([-xunit[1], xunit[0]])

   # print
   # print "Proplyd ", proplyd
   # print "D = ", D
   # print "xunit = ", xunit
   # print "yunit = ", yunit

x = []
y = []
R = []
theta = []
for shockpos in Shapes[proplyd]:
    # vector radius proplyd->shock normalized by separation D
    vecR = (shockpos - Centers[proplyd])/D
    R.append(np.hypot(*vecR))
    x.append(np.dot(vecR, xunit))
    y.append(np.dot(vecR, yunit))
    theta.append(np.arctan2(np.dot(vecR, yunit), np.dot(vecR, xunit)))

x = np.array(x)
y = np.array(y)
R = np.array(R)
theta = np.array(theta)


#Here is where the fit start
method_2 = "leastsq"
x_m = 0.0
y_m = 0.0
pixel_size = 0.05
# Assume 1-pixel uncertainty in positions
sigma = np.ones(len(x))*pixel_size


def circle_function(theta, xc, yc, rc):
    """Calculate distance from the origin as a function of theta for a
    circle that is centered at xc, yc, with radius rc

    """
    fac1 = xc*np.cos(theta) + yc*np.sin(theta)
    fac2 = xc*np.sin(theta) - yc*np.cos(theta)
    return fac1 + np.sqrt(rc**2 - fac2**2)


def circle_on_axis(theta, xc, rc):
    """Same as circle_function but with yc=0.
    """
    return circle_function(theta, xc, 0.0, rc)


def data_minus_model(thdata, rdata, rmodel_func, model_pars):
    return rdata - rmodel_func(thdata, *model_pars)

R_m = calc_R(x_m, y_m).mean()
if cmd_args.on_axis:
    # Two-parameter version
    p0 = x_m, R_m
    f = circle_on_axis
else:
    # Three-parameter version
    p0 = x_m, y_m, R_m
    f = circle_function

print "Initial parameter:"
print "Circle x, y, r = ", x_m, y_m, R_m
popt, pcov = optimize.curve_fit(f, theta, R, p0, sigma)

if cmd_args.on_axis:
    xc_2, R_2 = popt
    yc_2 = 0.0
else:
    xc_2, yc_2, R_2 = popt

print "Fit results:"
print "Circle x, y, r = ", xc_2, yc_2, R_2
print "Covariance matrix: ", pcov

# Deviation betwen model and data
error = data_minus_model(theta, R, f, popt)
print "Deviations between model and data:"
print "Max, rms: ", abs(error).max(), np.sqrt(np.mean(error**2))

theta_fit = np.linspace(-np.pi, np.pi, 180)
x_fit2 = xc_2 + R_2*np.cos(theta_fit)
y_fit2 = yc_2 + R_2*np.sin(theta_fit)

#***************************************************************************
#Plotting data
plt.plot(x, y, "o", label="{}: D = {:.2f} arcsec".format(proplyd, D))
#***************************************************************************
#Plotting the circle fit
plt.plot(x_fit2, y_fit2, 'k--', label=method_2, lw=2)
plt.plot([xc_2], [yc_2], 'gD', mec='r', mew=1)
#***************************************************************************

# Proplyd position (at the origin in this frame)
plt.plot(0.0, 0.0, "rx", label=None)

xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()

vmin = min(xmin, ymin)
vmax = max(xmax, ymax)

plt.xlim(xmin=vmin, xmax=vmax)
plt.ylim(ymin=vmin, ymax=vmax)

plt.legend(loc="best", prop=dict(size="x-small"))
plt.xlabel("z'/D'")
plt.ylabel("r'/R'")
plt.axis("equal")
plt.grid()
plt.title("{} fit circle".format(proplyd))
if cmd_args.on_axis:
    plt.savefig("LV-bowshocks-xy-onaxis-{}-{}.pdf".format(regfl_chsn,proplyd))
else:
    plt.savefig("LV-bowshocks-xy-{}-{}.pdf".format(regfl_chsn, proplyd))

