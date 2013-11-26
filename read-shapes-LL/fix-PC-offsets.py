import os
from astropy.io import fits
import json


def shift_fits_file(field,filters,dx,dy):
    """
    Shift all the PC images from an specific field
    """
    for f in filters:
        print "***** shifting {} image filter".format(f)
        print "dRA = {}, dDec = {}".format(dx,dy)
        hdulist = fits.open(os.path.join("../../HST",field+f+".fits"))
        hdulist[0].header["CRVAL1"]-=dx/3600.
        hdulist[0].header["CRVAL2"]-=dy/3600.
        outfile = "wcs_"+field+f+".fits"
        hdulist.writeto(os.path.join("../../HST",outfile))
        print "Newfile = {}".format(outfile)
    return


data = json.load(open("PC_offsets.json"))

# We can get all the PC fits files names by a combination of fields + filters + ".fits"
fields = ["1f","3f"]
filters = ["502","547","631","656","658","673"]

for field in fields:
    print "***** {} field *****"
    shift_fits_file(field,filters,data[field]["avdx"][0],data[field]["avdy"][0])
