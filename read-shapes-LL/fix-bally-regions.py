"""The final stage in correcting the Bally positions

Apply the corrections to the SOURCE-forma.reg and SOURCE-mask.reg
files

All the regions are of the form: 

SHAPE(RA,DEC[,OTHER_STUFF])

So we can just go through each line, correcting the RA and DEC
"""
import json
import glob
import os
import astropy.coordinates as coords
import astropy.units as u

data = json.load(open("bally-fieldshifts.json"))


def has_region(line):
    "True if line looks like a DS9 region"
    return "(" in line and ")" in line


def correct_coords(line, d, debug=False):
    """
    Apply corrections to RA and Dec in a ds9 region
    """
    ra, dec = line.split(",")[:2]
    ra = ra.split("(")[-1]
    dec = dec.split(")")[0]
    c = coords.ICRSCoordinates(ra, dec, unit=(u.hour, u.deg))
    new_ra = (c.ra - d["dx"]*u.arcsec).to_string(sep=":") 
    new_dec = (c.dec - d["dy"]*u.arcsec).to_string(sep=":")
    if debug:
        print "Replacing RA {} with {}".format(ra, new_ra)
        print "Replacing Dec {} with {}".format(dec, new_dec)
    return line.replace(ra, new_ra).replace(dec, new_dec)


for field, d in data.items(): 
    formas = glob.glob(os.path.join(field, "*-forma.reg"))
    masks = glob.glob(os.path.join(field, "*-mask.reg"))
    regfiles = formas + masks
    wcsfield = field.replace("_drz", "_wcs")
    if not os.path.isdir(wcsfield):
        os.mkdir(wcsfield)
    for regfile in regfiles:
        debug = False  # "LL1" in regfile
        newfile = regfile.replace("_drz", "_wcs")
        with open(regfile) as f:
            lines = f.readlines()
        newlines = []
        for line in lines:
            if has_region(line):
                newlines.append(correct_coords(line, d, debug))
            else:
                newlines.append(line)
        with open(newfile, "w") as f:
            f.writelines(newlines)
        
   
