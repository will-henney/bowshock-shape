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

def correct_coords(line, d):
    ra, dec = line.split(",")[:2]
    ra = ra.split("(")[-1]
    dec = ra.split(")")[0]
    dx = d["dx"]
    dy = d["dy"]
    c = coords.ICRSCoordinates(ra, dec, unit=(u.hour, u.deg))
    new_ra = (c.ra - dx*u.arcsec).to_string(sep=":") 
    new_dec = (c.dec - dy*u.arcsec).to_string(sep=":")
    return line.replace(ra, new_ra).replace(dec, new_dec)


for field, d in data.items(): 
    formas = glob.glob(os.path.join(field, "*-forma.reg"))
    masks = glob.glob(os.path.join(field, "*-mask.reg"))
    regfiles = formas + masks
    wcsfield = field.replace("_drz", "_wcs")
    if not os.path.isdir(wcsfield):
        os.mkdir(wcsfield)
    for regfile in regfiles:
        newfile = regfile.replace("_drz", "_wcs")
        with open(regfile) as f:
            lines = f.readlines()
        for line in lines:
            if has_region(line):
                line = correct_coords(line, d) + "\n"
        with open(newfile, "w") as f:
            f.writelines(lines)
        
   
