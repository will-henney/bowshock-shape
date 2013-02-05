import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(
    description=""" Choose a region file to work and the angle to measure radius""")

parser.add_argument("--region",type = str, 
                   choices = ("LV-OIII-positions-2.reg","LV-OIII-positions-3.reg","LV-Ha-positions-far.reg"),default = "LV-OIII-positions-3.reg",
                   help = " Choose a region file to work ")

parser.add_argument('--angle',type = float,default = 60, help = 'angle to measure radius')

cmd_args = parser.parse_args()
regionfile = cmd_args.region
angle = cmd_args.angle 
name,ext = regionfile.split('.')
regfl_chsn = name.split('-')[1] + name.split('-')[-1]

def extract_2nd(x):
    """
    This function extracts the second maximum value from an array, and also the reduced array, in case you want the third highest value or another.
    """
    xm = np.max(np.abs(x))
    xout =[]
    for y in x:
        if np.abs(y) != xm:
            xout.append(y)
#    return xout, np.max(xout)  #<----- This is a mistake. What'd happen if the max angle is negative?
    return xout,np.max(np.abs(xout))

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

print proplyds

for proplyd in proplyds: 
    # vector star->proplyd
    try:
        vecD = Centers['th1C'] - Centers[proplyd]
    except KeyError:
        # no center for this proplyd, so skip it
        continue
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
#        print "R, theta, x, y = {}, {}, {}, {}".format(R,theta, x, y)
#    th2nda,th2nd = extract_2nd(theta)
#    th3rda,th3rd = extract_2nd(th2nda)             #changing method
    theta_mean,rt = extract_n_mean(theta,R,angle)
    plt.plot(x, y, "o", label="{}: D = {:.2f} arcsec".format(proplyd, D))
    for x,y in zip(R,theta):
        if x == np.array(R).min():
            r0,th0 = x,y       
    print "{} {:.2f} {:.2f} {:.2f} {:.2f}".format(proplyd,r0,th0,rt,theta_mean)
    

plt.plot(0.0, 0.0, "rx", label=None)
plt.legend(loc="lower left", prop=dict(size="x-small"))
plt.xlabel("x")
plt.ylabel("y")
plt.axis("equal")
plt.grid()
plt.savefig("LV-bowshocks-xy-{}.pdf".format(regfl_chsn))

