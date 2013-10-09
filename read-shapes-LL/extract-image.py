import numpy as np
import json
from astropy.io import fits
from astropy import wcs
from astropy import coordinates as coord
from astropy import units as u 
import argparse

parser = argparse.ArgumentParser(
    description="""Extract a small image around a bowshock and save it as a FITS file""")

parser.add_argument("--source", type=str,
                    default="LL1",
                    help="Name of source (prefix for files)")

parser.add_argument("--fitsfile", type=str,
                    default="j8oc01010_drz.fits",
                    help="Original FITS image that contains the object")

parser.add_argument("--margin", type=float,
                    default=10.0,
                    help="Margin around object (arcseconds)")

parser.add_argument("--debug", action="store_true",
                    help="Print out verbose debugging info")

cmd_args = parser.parse_args()

hdu = fits.open(cmd_args.fitsfile)['SCI']
with open(cmd_args.source + "-xyc.json") as f:
     db = json.load(f)
w = wcs.WCS(hdu.header)

# w.wcs.print_contents()

star_pos = coord.ICRSCoordinates(
    ra=db["star"]["RA"], 
    dec=db["star"]["Dec"], 
    unit=(u.hour, u.degree)
)

##
## Construct a list of ra and dec of all points on either arc
##

def HMS(angle): 
    """
    Convert angle (which has astropy.units units) to an HMS string
    """
    return coord.Angle(angle).to_string(u.hour, sep=":")

def DMS(angle): 
    """
    Convert angle (which has astropy.units units) to a DMS string
    """
    return coord.Angle(angle).to_string(u.degree, sep=":")

ra = np.hstack([np.array(db[arc]["x"]) 
                for arc in "inner", "outer"]) * u.arcsecond + star_pos.ra
dec = np.hstack([np.array(db[arc]["y"]) 
                 for arc in "inner", "outer"]) * u.arcsecond + star_pos.dec

if cmd_args.debug:
    print HMS(ra)
    print DMS(dec)
##
## Find minimum and maximum RA, DEC
##

margin = cmd_args.margin * u.arcsecond
# Ignore the cos(delta) factor since it is ~= 1 for Orion
ra1 = ra.min() - margin
ra2 = ra.max() + margin
dec1 = dec.min() - margin
dec2 = dec.max() + margin

print "RA range: ", HMS(ra1), HMS(ra2)
print "Dec range: ", DMS(dec1), DMS(dec2)

##
## Rectangle in RA, Dec that encloses object with margin
##
coords = [
    [ra1.degree, dec1.degree], 
    [ra1.degree, dec2.degree], 
    [ra2.degree, dec1.degree], 
    [ra2.degree, dec2.degree],
]

##
## Convert to pixel coords and find enclosing rectangle
##
pix_coords = w.wcs_world2pix(coords, 0)
x = pix_coords[:,0]
y = pix_coords[:,1]
i1, i2 = int(x.min()), int(x.max()) + 1
j1, j2 = int(y.min()), int(y.max()) + 1

print "Extracted image window: [{}:{}, {}:{}]".format(i1, i2, j1, j2)

##
## Extract window from image and adjust WCS info
##
hdu.data = hdu.data[j1:j2, i1:i2]
hdu.header["CRPIX1"] -= i1
hdu.header["CRPIX2"] -= j1

##
## Save the small image
##
hdu.writeto(cmd_args.source + "-extract.fits", clobber=True)

