import numpy as np
import json
import os
from astropy.io import fits
from astropy import wcs
from astropy import coordinates as coord
from astropy import units as u 
import argparse
import fits_utils
import misc_utils

parser = argparse.ArgumentParser(
    description="""Extract a small image around a bowshock and save it as a FITS file""")

parser.add_argument("source", type=str,
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

image_name = os.path.basename(cmd_args.fitsfile).replace(".fits", "")
image_name = misc_utils.short_image_name(image_name)
hdu = fits_utils.get_image_hdu(fits.open(cmd_args.fitsfile), debug=cmd_args.debug)

dbfile = cmd_args.source + "-arcdata.json"
with open(dbfile) as f:
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
ra1 =  coord.Angle(ra.min() - margin ) 
ra2 =  coord.Angle(ra.max() + margin ) 
dec1 = coord.Angle(dec.min() - margin) 
dec2 = coord.Angle(dec.max() + margin) 

print "RA range: ", HMS(ra1), HMS(ra2)
print "Dec range: ", DMS(dec1), DMS(dec2)

##
## Rectangle in RA, Dec that encloses object with margin
##
coords = [
    [ra1.deg, dec1.deg], 
    [ra1.deg, dec2.deg], 
    [ra2.deg, dec1.deg], 
    [ra2.deg, dec2.deg],
]#.to(u.deg)


##
## Convert to pixel coords and find enclosing rectangle
##
pix_coords = w.wcs_world2pix(coords, 0)
x = pix_coords[:,0]
y = pix_coords[:,1]
i1, i2 = int(x.min()), int(x.max()) + 1
j1, j2 = int(y.min()), int(y.max()) + 1

ny, nx = hdu.data.shape
i1 = max(0, i1)
i2 = min(i2, nx-1)
j1 = max(0, j1)
j2 = min(j2, ny-1)
print "Extracted image window: [{}:{}, {}:{}]".format(i1, i2, j1, j2)

##
## Extract window from image and adjust WCS info
##
outhdu = fits.PrimaryHDU(
     data=hdu.data[j1:j2, i1:i2],
     header=hdu.header
)
outhdu.header["CRPIX1"] -= i1
outhdu.header["CRPIX2"] -= j1

# Only use header from astropy.wcs to correct
# non-standard stuff (e.g., equinox, dates)
whdr=w.to_header()
for kwd in "EQUINOX", "DATE-OBS":
     if kwd in whdr:
          if cmd_args.debug:
               print "Replacing {} = {} with {}".format(
                    kwd, outhdu.header[kwd], whdr[kwd]
               )
          outhdu.header[kwd] = whdr[kwd]

## Further correction to EQUINOX keyword if necessary
equinox = outhdu.header.get("EQUINOX")
if isinstance(equinox, basestring):
     if cmd_args.debug:
          print "Converting EQUINOX to float"
     outhdu.header["EQUINOX"] = 2000.0

## TODO: copy over more keywords from the original FITS header
## E.g., filter, camera, etc

##
## Save the small image
##
outfile = "-".join([cmd_args.source, image_name, "extract.fits"])
outhdu.writeto(outfile, output_verify="fix", clobber=True)


##
## Update the JSON database with information about the image
##

# Always overwrite an existing section with the same name
db[image_name] = {
    "original FITS file": cmd_args.fitsfile,
    "extracted FITS file": outfile,
}
db[image_name].update(fits_utils.get_instrument_configuration(hdu))
db["info"]["history"].append("Image " + image_name + " added by " + misc_utils.run_info())

misc_utils.update_json_file(db, dbfile)
