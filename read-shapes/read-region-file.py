import sys
import numpy as np
import matplotlib.pyplot as plt
regionfile = "LV-OIII-positions-2.reg" 


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

    print
    print "Proplyd ", proplyd
    print "D = ", D
    print "xunit = ", xunit
    print "yunit = ", yunit

    x = []; y = []; R = []; theta = []
    for shockpos in Shapes[proplyd]:
        # vector radius proplyd->shock normalized by separation D
        vecR = (shockpos - Centers[proplyd])/D
        R.append(np.hypot(*vecR))
        x.append(np.dot(vecR, xunit))
        y.append(np.dot(vecR, yunit))
        theta.append(np.degrees(np.arctan2(y, x)))
        # print "R, theta, x, y = {}, {}, {}, {}".format(R, np.degrees(theta), x, y)

    plt.plot(x, y, "o", label="{}: D = {:.2f} arcsec".format(proplyd, D))

plt.plot(0.0, 0.0, "x", label=None)
plt.legend(loc="lower left")
plt.xlabel("x")
plt.ylabel("y")
plt.axis("equal")
plt.grid()
plt.savefig("LV-bowshocks-xy.pdf")

