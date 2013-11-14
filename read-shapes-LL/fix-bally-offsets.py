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
    "I": ["01", "16", "14", "17", "06", "24"],
    "II": ["07", "08"],
    "III": ["09"],
    "IV": ["02"],
    "V": ["18"],
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
        files.append(data[field]["file"] + ".fits")
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
        hdu = fits.open(os.path.join(fits_dir, fn))["SCI"]
    except IOError:
        print fn, "not found"
        return
    hdu.header["CRVAL1"] += dx/3600.0
    hdu.header["CRVAL2"] += dy/3600.0
    outfile = fn.replace(*replace)
    hdu.writeto(os.path.join(fits_dir, outfile), clobber=True)
    return

# Now actually update the FITS files
for data in groupdata.values():
    for fn in data["files"]:
        shift_fits_file(fn, data["dx"], data["dy"])
