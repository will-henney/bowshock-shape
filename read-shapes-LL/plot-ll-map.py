from __future__ import print_function
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import json
import os


nicknames = {
    "166-316": "LV2b", 
    "177-341": "HST1", 
    "167-317": "LV2",  
    "163-317": "LV3",  
    "161-324": "LV4",  
    "158-323": "LV5", 
}

def find(name, path):
    """
    Original from http://stackoverflow.com/questions/1724693/find-a-file-in-python
    """
    for root, dirs, files in os.walk(path):
        for realname in name, "w" + name:
            if realname in files and not "_drz" in root:
                return os.path.join(root, realname)
    return None
            

def plot_map(limits, figname, canvas_size, innerbox=None, arrowscale=1.0):
    plt.clf()
    c = coord.SkyCoord(RAs, Decs, unit=(u.hourangle, u.degree))
    x, y = c.ra.deg, c.dec.deg
    x = -(x - x0)*np.cos(c.dec.radian)*3600.0
    y = (y - y0)*3600.0

    print(x, y)
    plt.plot(x, y, "o", alpha=0.2)
    for label, xx, yy in zip(names, x, y):
        
        #
        # Try and draw the inner and outer arcs
        #

        name = label.split()[-1]
        if name in problem_sources:
            continue

        nickname = nicknames.get(name, name)
        jsonfile = "{}-arcdata.json".format(name)
        jsonfilex = "{}-arcdata.json".format(nickname)
        found = find(jsonfilex, "../JorgeBowshocks/Jorge_prop/PC-will")
        if found is not None:
            f = open(found)
        else:
            found = find(jsonfile, ".")
            if found is not None:
                f = open(found)
            else:
                print("Could not open ", jsonfile)
                continue
        arc_data = json.load(f)

        # Second, load in the data and draw the arcs
        small_A = []
        for arc, color in ["inner", "m"], ["outer", "g"]:
            if arc in arc_data:
                dx = np.array(arc_data[arc]["x"])
                dy = np.array(arc_data[arc]["y"])
                plt.plot(xx - dx, yy + dy, "-" + color, lw=1.0, alpha=0.6)
                print("Plotted {} arc for {}".format(arc, found))
                if "Rc" in arc_data[arc]:
                    xc = arc_data[arc]["xc"]
                    yc = arc_data[arc]["yc"]
                    Rc = arc_data[arc]["Rc"]
                    PAc = np.radians(arc_data[arc]["PAc"])
                    if arc_data[arc]["Rc"] < arc_data[arc]["R0"]:
                        # Flip the arrow for Rc < R0
                        PAc += np.pi
                    if arc_data[arc]["Rc"] < 1.5*arc_data[arc]["R0"]:
                        small_A.append(arc)
                        PAm = np.radians(arc_data[arc]["PA0"]
                                         + np.mean(arc_data[arc]["theta"]))
                    else:
                        PAm = None
                    # Plot the fitted circle if present
                    plt.plot(xx - xc, yy + yc, "+k", ms=2.0)
                    c = plt.Circle((xx - xc, yy + yc), radius=Rc, fc='none', ec="k", alpha=0.2, lw=0.2)
                    plt.gca().add_patch(c)
                    PA = PAc if PAm is None else PAm
                    ax = -0.5*Rc*np.sin(PA)
                    ay = 0.5*Rc*np.cos(PA)
                    plt.arrow(xx-xc, yy+yc, 4*ax*arrowscale, 4*ay*arrowscale,
                              fc='none', ec=color, 
                              width=0.001, alpha=0.8, lw=1.5,
                              head_width=2.0*arrowscale, head_length=4.0*arrowscale,
                          )

        if innerbox is None:
            skip_annotation = False
        else:
            x1, x2, y1, y2 = innerbox
            skip_annotation = (x1 <= xx <= x2) and (y1 <= yy <= y2)
        
        # Order of octants is anticlockwise around square starting at top:
        #    1 0 7
        #    2 * 6
        #    3 4 5
        alignment_by_octant = [
            ('center', 'bottom'),
            ('right', 'bottom'),
            ('right', 'center'),
            ('right', 'top'),
            ('center', 'top'),
            ('left', 'top'),
            ('left', 'center'),
            ('left', 'bottom'),
        ]
        if not skip_annotation:
            boxcolor = 'orange' if name in interprop_sources else 'white'
            if len(small_A) == 2:
                labelcolor = 'red'
            elif 'inner' in small_A:
                labelcolor = 'magenta'
            elif 'outer' in small_A:
                labelcolor = 'green'
            else:
                labelcolor = 'black'
            PA = np.radians(arc_data["star"]["PA"] + 180.0)
            ioctant = int(((np.degrees(PA) + 22.5) % 360)*8/360)
            print('Octant:', ioctant, 'PA:', np.degrees(PA))
            ha, va = alignment_by_octant[ioctant]
            xytext = (-3*np.sin(PA), 3*np.cos(PA))
            plt.annotate(label, (xx, yy), alpha=0.8, size=5, color=labelcolor,
                         xytext=xytext, textcoords='offset points',
                         ha=ha, va=va,
                         bbox={'facecolor': boxcolor, 
                               'alpha': 0.5,
                               'pad': 2,
                               'linewidth': 0.1,
                           },
            )

    c = coord.SkyCoord(pRAs, pDecs, unit=(u.hourangle, u.degree))
    x, y = c.ra.deg, c.dec.deg
    x = -(x - x0)*np.cos(c.dec.radian)*3600.0
    y = (y - y0)*3600.0
    plt.scatter(x, y, c=pColors, s=pSizes, edgecolors='none', alpha=0.5, zorder=100)

    if innerbox is not None:
        plt.fill_between([x1, x2], [y1, y1], [y2, y2],
                         edgecolor='none', facecolor='k', alpha=0.1)
    plt.axis("equal")
    plt.axis(limits)
    plt.xlabel("RA offset, arcsec")
    plt.ylabel("Dec offset, arcsec")
    plt.title("Positions of LL objects in the Orion Nebula")
    plt.grid()
    plt.gcf().set_size_inches(canvas_size)
    plt.savefig(figname)


if __name__ == "__main__":

    #
    # Set up arc data
    #

    table = Table.read("luis-programas/arcs-summary-merge.tab", 
                     format="ascii.commented_header", delimiter="\t",
                     fill_values=('--', np.nan) ).filled(np.nan)
    names = table["Object"].data
    RAs = table["RA"].data
    Decs = table["Dec"].data

    # TABLE_FILE = "ll-data.json"
    # table = json.load(open(TABLE_FILE))

    # names = table.keys()
    # RAs = [v["RA"] for v in table.values()]
    # Decs = [v["Dec"] for v in table.values()]
    # Fields = [v["Bally"] for v in table.values()]

    #
    # Set up proplyd data 
    #
    pColor_from_Type = {
        "i": "r", "d": "k", "rn": "c", "j": "g",
    }
    pSize_from_Type = {
        "i": 1.0, "d": 2.0, "rn": 0.7, "j": 0.7,
    }

    PROPLYD_TABLE_FILE = "ricci-data.json"
    ptable = json.load(open(PROPLYD_TABLE_FILE))

    pRAs = [v["RA"] for v in ptable.values()]
    pDecs = [v["Dec"] for v in ptable.values()]
    pColors = [pColor_from_Type[v["Type"]] for v in ptable.values()]
    pSizes = [pSize_from_Type[v["Type"]] for v in ptable.values()]
    pSizes = np.array(pSizes)*15.0

    with open("luis-programas/problem-sources.txt") as f:
        problem_sources = f.read().split('\n')
    with open("luis-programas/interproplyd.txt") as f:
        interprop_sources = f.read().split('\n')

    c0 = coord.SkyCoord("05:35:16.463", "-05:23:23.18",
                        unit=(u.hourangle, u.degree))
    x0, y0 = c0.ra.deg, c0.dec.deg

    zoombox = [-50, 50, -35, 65]
    fullbox = [-350, 600, -650, 200]
    plot_map(fullbox, "ll-positions.pdf", (10, 10), innerbox=zoombox, arrowscale=2.0)
    plot_map(zoombox, "ll-positions-zoom.pdf", (10, 10), arrowscale=0.7)










