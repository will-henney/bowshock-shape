from __future__ import print_function
import json
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.units as u
import astropy.coordinates as coord
import aplpy
from misc_utils import expand_fits_path

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
            

def plot_map(limits, figname, canvas_size,
             fitsfile='$LARGE_FITS_DIR/wfi/Orion_H_A_shallow.fits',
             north=False,
             vmin=0.0, vmax=None, stretch='linear',
             innerbox=None, arrowscale=1.0):
    # Use an image as a backdrop
    fig = aplpy.FITSFigure(expand_fits_path(fitsfile),
                           figsize=canvas_size, north=north)
    # Set the viewport
    xc, yc = (limits[0] + limits[1])/2, (limits[2] + limits[3])/2
    w, h = limits[1] - limits[0], limits[3] - limits[2]
    fig.recenter(c0.ra.deg - xc/3600, c0.dec.deg + yc/3600,
                 width=w/3600, height=h/3600)
    fig.show_grayscale(vmin=vmin, vmax=vmax, invert=True,
                       stretch=stretch, interpolation='none')
    ax = fig._ax1


    c = coord.SkyCoord(RAs, Decs, unit=(u.hourangle, u.degree))
    # Cartesian pixel coordinates of each source
    x, y = fig.world2pixel(c.ra.deg, c.dec.deg)

    # Pixel size in degrees
    pix_scale = aplpy.wcs_util.celestial_pixel_scale(fig._wcs)
    # Convert to arcsec
    pix_scale *= 3600

    ax.plot(x, y, "o", alpha=0.2)
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
        for arc, color in ["inner", "m"], ["outer", "g"]:
            if arc in arc_data:
                dx = np.array(arc_data[arc]["x"])/pix_scale
                dy = np.array(arc_data[arc]["y"])/pix_scale
                ax.plot(xx - dx, yy + dy, "-" + color, lw=1.0, alpha=0.6)
                print("Plotted {} arc for {}".format(arc, found))
                if "Rc" in arc_data[arc]:
                    xc = arc_data[arc]["xc"]/pix_scale
                    yc = arc_data[arc]["yc"]/pix_scale
                    Rc = arc_data[arc]["Rc"]/pix_scale
                    PAc = np.radians(arc_data[arc]["PAc"])
                    PAm = np.radians(arc_data[arc]["PA0"]
                                     + np.mean(arc_data[arc]["theta"]))

                    # Plot the fitted circle if present
                    ax.plot(xx - xc, yy + yc, "+k", ms=2.0)
                    c = plt.Circle((xx - xc, yy + yc), radius=Rc, fc='none', ec="k", alpha=0.2, lw=0.2)
                    ax.add_patch(c)
                    PA = PAm if arc_data[arc]["Rc"] < 1.5*arc_data[arc]["R0"] else PAc
                    arrowx = -0.5*Rc*np.sin(PA)
                    arrowy = 0.5*Rc*np.cos(PA)
                    ax.arrow(xx-xc, yy+yc, 4*arrowx*arrowscale, 4*arrowy*arrowscale,
                             fc='none', ec=color, 
                             width=0.001, alpha=0.8, lw=1.5,
                             head_width=8.0*arrowscale, head_length=16.0*arrowscale,
                         )

        if innerbox is None:
            skip_annotation = False
        else:
            x1, y1 = fig.world2pixel(c0.ra.deg - innerbox[0]/3600,
                                     c0.dec.deg + innerbox[2]/3600)
            x2, y2 = fig.world2pixel(c0.ra.deg - innerbox[1]/3600,
                                     c0.dec.deg + innerbox[3]/3600)
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
            PA = np.radians(arc_data["star"]["PA"] + 180.0)
            ioctant = int(((np.degrees(PA) + 22.5) % 360)*8/360)
            print('Octant:', ioctant, 'PA:', np.degrees(PA))
            ha, va = alignment_by_octant[ioctant]
            xytext = (-3*np.sin(PA), 3*np.cos(PA))
            ax.annotate(label, (xx, yy), alpha=0.8, size=5,
                        xytext=xytext, textcoords='offset points',
                        ha=ha, va=va,
                        bbox={'facecolor': boxcolor, 
                              'alpha': 0.5,
                              'pad': 2,
                              'linewidth': 0.1,
                          },
            )

    c = coord.SkyCoord(pRAs, pDecs, unit=(u.hourangle, u.degree))
    x, y = fig.world2pixel(c.ra.deg, c.dec.deg)
    ax.scatter(x, y, c=pColors, s=pSizes, edgecolors='none', alpha=0.5, zorder=100)

    if innerbox is not None:
        ax.fill_between([x1, x2], [y1, y1], [y2, y2],
                         edgecolor='black', lw=0.5, facecolor='yellow', alpha=0.3)

    fig.save(figname)


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

    zoombox = [-50, 50, -35, 65]
    fullbox = [-350, 600, -700, 250]
    plot_map(fullbox, "ll-pos-image.pdf", (10, 10),
             vmin=5.0, vmax=2000.0, stretch='sqrt',
             innerbox=zoombox, arrowscale=2.0)
    plot_map(zoombox, "ll-pos-image-zoom.pdf", (10, 10),
             fitsfile='$LARGE_FITS_DIR/acs/hlsp_orion_hst_acs_strip0l_f658n_v1_drz.fits',
             north=False,
             vmin=8.0, vmax=400.0,
             arrowscale=0.7)










