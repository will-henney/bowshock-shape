import json
import sys
import os
from astropy.io import fits
import image_statistics

try:
    fits_dir = sys.argv[1]
except: 
    sys.exit("Usage: {} FITS_DIR".format(sys.argv[0]))

# Groups of fields that seem to have common offsets
groups = {
    "I": ["01", "16", "17", "06", "24", "07"],
    "II": ["08"],
    "III": ["09"],
    "IV": ["02"],
    "V": ["18"],
    "VI": ["14"]
}

data = json.load(open("bally-offsets.json"))

# Aggregate offsets for each group and find average
groupdata = {}
for group, fields in groups.items():
    dx = []
    dy = []
    files = []
    for field in fields:
        dx.extend(data[field]["dx"])
        dy.extend(data[field]["dy"])
        files.append(data[field]["file"])
    avdx, sigdx = image_statistics.trimean_and_iqr(dx)
    avdy, sigdy = image_statistics.trimean_and_iqr(dy)
    groupdata[group] = {
        "files": files,
        "dx": avdx,
        "dy": avdy,
        "N": len(dx)
    }

def shift_fits_file(fn, dx, dy, replace=("_drz", "_wcs")):
    try:
        hdulist = fits.open(os.path.join(fits_dir, fn))
    except IOError:
        print fn, "not found"
        return
    hdulist["SCI"].header["CRVAL1"] -= dx/3600.0
    hdulist["SCI"].header["CRVAL2"] -= dy/3600.0
    outfile = fn.replace(*replace)
    hdulist.writeto(os.path.join(fits_dir, outfile), clobber=True)
    return

# Now actually update the FITS files
fieldshifts = {}
for data in groupdata.values():
    for fn in data["files"]:
        shift_fits_file(fn + ".fits", data["dx"], data["dy"])
        fieldshifts[fn] = {"dx": data["dx"], "dy":  data["dy"]}

# And save a dict of the shifts for each field
# so we ca use them later to correct the .reg files
with open("bally-fieldshifts.json", "w") as f:
    json.dump(fieldshifts, f, indent=2)
