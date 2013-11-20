from astropysics import coords
import numpy as np
import matplotlib.pyplot as plt
import json
import os

def find(name, path):
    """
    Original from http://stackoverflow.com/questions/1724693/find-a-file-in-python
    """
    for root, dirs, files in os.walk(path):
        if name in files and not "_drz" in root:
            return os.path.join(root, name)
    return None
            

def plot_map(limits, figname, canvas_size, exclude=0.0):
    x, y = coords.radec_str_to_decimal(RAs, Decs) 
    x = -(x - x0)*3600.0
    y = (y - y0)*3600.0

    plt.plot(x, y, "o", alpha=0.2)
    for label, xx, yy, field_list in zip(names, x, y, Fields):
        
        if xx*2 + yy**2 > exclude**2:
            plt.annotate(label, (xx, yy), alpha=0.8, size=5,
                         xytext=(-2,2), textcoords='offset points',
                         ha='right', va='bottom',
                         bbox={'facecolor': 'white', 
                               'alpha': 0.5,
                               'pad': 2,
                               'linewidth': 0.1,
                           },
            )
        #
        # Try and draw the inner and outer arcs
        #

        # First, work out where the JSON file is...
        try:
            field = "{:02d}".format(field_list)
        except:
            try:
                field = field_list[:2]
            except:
                print "failed"
                continue

        # folder = "j8oc{}010_drz".format(field)
        jsonfile = "{}-arcdata.json".format(label.split()[-1])
        found = find(jsonfile, "..")
        if found is not None:
            f = open(found)
        else:
            print "Could not open ", jsonfile
            continue

        # Second, load in the data and draw the arcs
        arc_data = json.load(f)
        for arc, color in ["inner", "m"], ["outer", "g"]:
            if arc in arc_data:
                dx = np.array(arc_data[arc]["x"])
                dy = np.array(arc_data[arc]["y"])
                plt.plot(xx - dx, yy + dy, "-" + color, lw=1.0, alpha=0.6)
                print "Plotted {} arc for {}".format(arc, found)
                # Now try and draw arrows for the arc axes too
                R0 = arc_data[arc]["R0"]
                PA = np.radians(arc_data[arc]["PA0"])
                ax = -R0*np.sin(PA)
                ay = R0*np.cos(PA)
                plt.arrow(xx, yy, 4*ax, 4*ay, fc='none', ec=color, 
                          width=0.0003, alpha=0.6, lw=0.5,
                          head_width=2.0, head_length=4.0,
                )
                if "Rc" in arc_data[arc]:
                    xc = arc_data[arc]["xc"]
                    yc = arc_data[arc]["yc"]
                    Rc = arc_data[arc]["Rc"]
                    PAc = np.radians(arc_data[arc]["PAc"])
                    # Plot the fitted circle if present
                    plt.plot(xx - xc, yy + yc, "+k", ms=2.0)
                    c = plt.Circle((xx - xc, yy + yc), radius=Rc, fc='none', ec="k", alpha=0.2, lw=0.2)
                    plt.gca().add_patch(c)
                    ax = -0.5*Rc*np.sin(PAc)
                    ay = 0.5*Rc*np.cos(PAc)
                    plt.arrow(xx-xc, yy+yc, 4*ax, 4*ay, fc='none', ec=color, 
                              width=0.0003, alpha=0.3, lw=0.3,
                              head_width=2.0, head_length=4.0,
                          )


    x, y = coords.radec_str_to_decimal(pRAs, pDecs)
    x = -(x - x0)*3600.0
    y = (y - y0)*3600.0
    plt.scatter(x, y, c=pColors, s=pSizes, faceted=False, alpha=0.5, zorder=100)

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
    TABLE_FILE = "ll-data.json"
    table = json.load(open(TABLE_FILE))

    names = table.keys()
    RAs = [v["RA"] for v in table.values()]
    Decs = [v["Dec"] for v in table.values()]
    Fields = [v["Bally"] for v in table.values()]


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


    [x0], [y0] = coords.radec_str_to_decimal(["05:35:16.463"], ["-05:23:23.18"])

    plot_map([-350, 600, -650, 200], "ll-positions.pdf", (20, 20), exclude=50.0)
    plot_map([-50, 50, -35, 65], "ll-positions-zoom.pdf", (10, 10))










