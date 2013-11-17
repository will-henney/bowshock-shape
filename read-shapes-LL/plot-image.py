import json
from astropy.io import fits
from astropy import coordinates as coord
from astropy import units as u 
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import aplpy
import pyregion
from misc_utils import run_info, update_json_file
from fits_utils import get_image_hdu

parser = argparse.ArgumentParser(
    description="""Plot side-by-side pure image of source and image with overlays""")

parser.add_argument("source", type=str,
                    help="Name of source (prefix for files)")
parser.add_argument("--image", type=str,
                    default="j8oc01010_drz",
                    help="Name of original FITS image (section in database)")
parser.add_argument("--maxfactor", type=float, default=None,
                    help="""\
Set the maximum brightness in units of shell dispersions above shell average.
Default is to use any value already saved in the JSON file
with fall-back of 3.0""")
parser.add_argument("--minfactor", type=float, default=None,
                    help="""\
Set the minimum brightness in units of shell dispersions below shell average.
Default is to use any value already saved in the JSON file
with fall-back of 3.0""")
parser.add_argument("--debug", action="store_true",
                    help="Print out verbose debugging info")

MAXFACTOR_DEFAULT = 3.0
MINFACTOR_DEFAULT = 3.0
cmd_args = parser.parse_args()

dbfile = cmd_args.source + "-arcdata.json"
arcdata = json.load(open(dbfile))

image_name = cmd_args.image

if not image_name in arcdata:
    raise ValueError, image_name + " not found - try running extract-image.py first"
if not "shell" in arcdata[image_name]:
    raise ValueError, "Shell data not found - try running arc_brightness.py first"



fitsfile = arcdata[image_name]["extracted FITS file"]
hdulist = fits.open(fitsfile)
hdu = get_image_hdu(hdulist, debug=cmd_args.debug)


# Find coordinates of the central star and radius of curvature
# We want all the values in degrees for use with aplpy
ra0 = coord.Longitude(arcdata["star"]["RA"], unit=u.hour).to(u.deg).value
dec0 = coord.Latitude(arcdata["star"]["Dec"], unit=u.deg).value
Rc = arcdata["outer"]["Rc"] * u.arcsec / u.deg


#
# Find brightness limits for plot
#

# Average and dispersion for shell and background
avsh = max(arcdata[image_name]["shell"]["value"], 
           arcdata[image_name]["shell center"]["value"])
dsh = max(arcdata[image_name]["shell"]["delta"], 
          arcdata[image_name]["shell center"]["delta"])
avbg = arcdata[image_name]["background"]["value"]
dbg = arcdata[image_name]["background"]["delta"]

plot_limits_modified = False

# Maximum
vmax = None
if cmd_args.maxfactor is None:
    if "plot limits" in arcdata[image_name]:
        # Use previously stored value if present
        vmax = arcdata[image_name]["plot limits"]["max"]
    else:
        maxfactor = MAXFACTOR_DEFAULT
else:
    maxfactor = cmd_args.maxfactor
if vmax is None:
    vmax = avsh + maxfactor*dsh
    plot_limits_modified = True
    
# Minimum
vmin = None
if cmd_args.minfactor is None:
    if "plot limits" in arcdata[image_name]:
        # Use previously stored value if present
        vmin = arcdata[image_name]["plot limits"]["min"]
    else:
        minfactor = MINFACTOR_DEFAULT
else:
    minfactor = cmd_args.minfactor
if vmin is None:
    vmin = avbg - minfactor*dbg
    plot_limits_modified = True

# Now update the db if necessary
if plot_limits_modified:
    arcdata[image_name]["plot limits"] = {"max": vmax, "min": vmin}
    arcdata["info"]["history"].append(
        "Plot limits modified for " + image_name + " " + run_info())
    arcdata["help"]["brightness"]["plot limits"] \
        = "Minimum and maximum brightness for *-images.pdf"
    update_json_file(arcdata, dbfile)
    
#
# Plot image of the FITS array of this object
# 
plt.clf()
f = plt.figure(figsize=(18,9))

ax1 = aplpy.FITSFigure(hdu, figure=f, subplot=(1, 2, 1), north=True)
ax1.recenter(ra0, dec0, 2*Rc)
ax1.show_grayscale(invert=True, vmin=vmin, vmax=vmax)
ax1.add_colorbar()
ax1.colorbar.hide()             # necessary to get panels to be same size

ax2 = aplpy.FITSFigure(hdu, figure=f, subplot=(1, 2, 2), north=True)
ax2.recenter(ra0, dec0, 2*Rc)
ax2.show_grayscale(invert=True, vmin=vmin, vmax=vmax)
ax2.show_regions(cmd_args.source + "-forma.reg")
ax2.show_regions(cmd_args.source + "-arcfits.reg")
# With the mask regions, the file may not exist
try:
    mask_regions = pyregion.open(cmd_args.source + "-mask.reg")
except IOError:
    print "No mask regions found"
else:
    # If it does, then we switch the color to yellow and linestyle to dashed
    for r in mask_regions:
        r.attr[1]['color'] = 'yellow'
        r.attr[1]['dash'] = '1 '
    ax2.show_regions(mask_regions)
ax2.add_scalebar(5.0/3600)
ax2.scalebar.set_label('5 arcsec')
ax2.scalebar.set(color='orange', linewidth=4, alpha=0.9)
ax2.scalebar.set_font(size='medium', weight='bold',
                      stretch='normal', family='sans-serif',
                      style='normal', variant='normal')
ax2.add_label(0.5, 0.3, cmd_args.source, color="white",
              weight='bold', size='x-large', relative=True, zorder=1000)
dx, dy = 0.001, -0.001
ax2.add_label(0.5+dx, 0.3+dy, cmd_args.source, color="black", alpha=0.6,
              weight='bold', size='x-large', relative=True, zorder=999)
ax2.axis_labels.hide_y()
ax2.tick_labels.hide_y()
ax2.add_colorbar()
ax2.colorbar.set_axis_label_text("counts")


f.tight_layout()
f.savefig("-".join([cmd_args.source, image_name, "images.pdf"]))

# record the --maxfactor and the --minfactor in the *-xycb.json file
# also their respective help section











