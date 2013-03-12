import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy import optimize

parser = argparse.ArgumentParser(
    description=""" Choose a region file to work and the angle to measure radius""")

parser.add_argument("--region",type = str, 
                   choices = ("LV-OIII-positions-2.reg","LV-OIII-positions-3.reg","LV-Ha-positions-far.reg"),default = "LV-OIII-positions-3.reg",
                   help = " Choose a region file to work ")

#parser.add_argument('--angle',type = float,default = 60, help = 'angle to measure radius')
parser.add_argument('--proplyd',type = str, default= 'LV3',help = 'choose a proplyd to work')
cmd_args = parser.parse_args()
regionfile = cmd_args.region
#angle = cmd_args.angle 
name,ext = regionfile.split('.')
regfl_chsn = name.split('-')[1] + name.split('-')[-1]

def calc_R(xc, yc):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((x-xc)**2 + (y-yc)**2)

#@countcalls
def f_2(c):
    """ calculate the algebraic distance between the 2D points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(*c)
    return Ri - Ri.mean()

def wf_2(c):
    """ calculate the algebraic distance between the 2D weighted points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(*c)
    return u*(Ri - Ri.mean())

def extract_n_mean(A,B,value):
    """
    This function extracts the elements from an array with numerical value near to the input "value", and gets the mean of those elements. If the ith
    element is extracted from array A, will do the same with array B
    """
    interval = []
    intervalr = []
    for x,y in zip(A,B):
        if ( (np.abs(x) <= value +5) & (np.abs(x) >= value-5) ):
            interval.append(np.abs(x))
            intervalr.append(y)
    exAmean = np.mean(interval)
    exBmean = np.mean(intervalr)
    return exAmean,exBmean

def extract_data(line):
    """
    Find the coordinates and description text from the line of a ds9 regions file
    """
    coordstring, paramstring = line.split("#")
    shape, numbers = coordstring.split(")")[0].split("(")
    ra, dec = numbers.split(",")[0:2]
    if "tag" in line:
        text = paramstring.split("tag={")[1].split("}")[0]
    elif "text" in line:
        text = paramstring.split("text={")[1].split("}")[0]
    else:
        text = "NO TEXT"
    return ra, dec, text

Shapes = {}
Centers = {}

with open(regionfile) as f:
    lines = f.readlines()
    for line in lines: 
        skipthisline = line.startswith("#") \
            or line.startswith("global") \
            or not "#" in line 
        if skipthisline: 
            continue
        ra, dec, text = extract_data(line)
        # name of source (LV knot or star)
        source = text.split()[0]
        if "OW" in text or text == "NO TEXT":
            # We do not use the OW positions
            continue

        shr,smin,ssec = ra.split(':')
        hr,mn,sec = float(shr),float(smin),float(ssec)
        ra_arcsec = 15*(3600.*hr+60.*mn+sec)
        sdeg,samin,sasec = dec.split(':')
        deg,amin,asec = float(sdeg),float(samin),float(sasec)
        dec_arcsec = 3600.*deg + np.sign(deg)*(60*amin + asec)

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

x = []; y = []; R = []; theta = []
for shockpos in Shapes[proplyd]:
    # vector radius proplyd->shock normalized by separation D
    vecR = (shockpos - Centers[proplyd])/D
    R.append(np.hypot(*vecR)) 
    x.append(np.dot(vecR, xunit))
    y.append(np.dot(vecR, yunit))
    theta.append( np.degrees( np.arctan2( np.dot(vecR,yunit),np.dot(vecR,xunit) ) ) )
#    print "R, theta, x, y = {}, {}, {}, {}".format(R,theta, x, y)

#**********************************************************************************************

#Here is where the fit start
method_2  = "leastsq"
x_m = np.mean(x)
y_m = np.mean(y)
u = np.ones(len(x))


center_estimate = x_m, y_m
center_2, ier = optimize.leastsq(wf_2, center_estimate)


xc_2, yc_2 = center_2
Ri_2       = calc_R(xc_2, yc_2)
R_2        = Ri_2.mean()
residu_2   = sum((Ri_2 - R_2)**2)
residu2_2  = sum((Ri_2**2-R_2**2)**2)
theta_fit = np.linspace(-np.pi, np.pi, 180)
x_fit2 = xc_2 + R_2*np.cos(theta_fit)
y_fit2 = yc_2 + R_2*np.sin(theta_fit)

#**********************************************************************************************
#Calculate R_45 and R_0
#theta_mean,rt = extract_n_mean(theta,R,45)
for r,t in zip(R,theta):
    if r == np.array(R).min():
        r0,th0 = r,t


#**********************************************************************************************
#Plotting data
plt.plot(x, y, "o", label="{}: D = {:.2f} arcsec".format(proplyd, D))
#**********************************************************************************************
#Plotting R_45 (both sides)
#plt.plot(rt*np.cos(np.radians(theta_mean)),rt*np.sin(np.radians(theta_mean)),'*',label = None)
#plt.plot(rt*np.cos(np.radians(theta_mean)),-rt*np.sin(np.radians(theta_mean)),'*',label = None)
#**********************************************************************************************
#Calculating R_60, plotting R_60 and R_0
#theta_mean,rt = extract_n_mean(theta,R,60)
#plt.plot(rt*np.cos(np.radians(theta_mean)),rt*np.sin(np.radians(theta_mean)),'p',label = None)
#plt.plot(rt*np.cos(np.radians(theta_mean)),-rt*np.sin(np.radians(theta_mean)),'p',label = None)
#plt.plot(r0*np.cos(np.radians(th0)),r0*np.sin(np.radians(th0)),'s',label = None)
#**********************************************************************************************
#Plotting the circle fit
plt.plot(x_fit2, y_fit2, 'k--', label=method_2, lw=2)
plt.plot([xc_2], [yc_2], 'gD', mec='r', mew=1)
#**********************************************************************************************
print "{} {:.2f} {:.2f} {:.2f}".format(proplyd,r0,th0,R_2)
    

plt.plot(0.0, 0.0, "rx", label=None) # Proplyd position (at the origin in this frame)

xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()

vmin = min(xmin, ymin)
vmax = max(xmax, ymax)

plt.xlim(xmin=vmin, xmax=vmax)
plt.ylim(ymin=vmin, ymax=vmax)

plt.legend(loc="best", prop=dict(size="x-small"))
plt.xlabel("x")
plt.ylabel("y")
plt.axis("equal")
plt.grid()
plt.title("{} fit circle".format(proplyd))
plt.savefig("LV-bowshocks-xy-{}-{}.pdf".format(regfl_chsn,proplyd))

