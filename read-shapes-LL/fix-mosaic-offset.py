import os
import json
import numpy as np
import sys
from astropy.io import fits
import image_statistics
import glob

fits_dir = "../HST"

data = json.load(open("mosaic_offsets.json"))

mosaic_pattern = "GO5469PCf*.fits"
path_fits = os.path.join(fits_dir,mosaic_pattern)

#print glob.glob(path_fits) #This function is just what I needed

def shift_fits_file(dx,dy,fts):        
    hdulist = fits.open(fts)
    hdulist[0].header["CRVAL1"] -=dx/3600.
    hdulist[0].header["CRVAL2"] -=dy/3600.
    outfile = fts.replace("GO","wcs_GO")
    hdulist.writeto(outfile)
    return

for fitsfile in glob.glob(path_fits):
    print "*****Shifting {} fits file*****".format(fitsfile)
    shift_fits_file(data["avdx"][0],data["avdy"][0],fitsfile)
