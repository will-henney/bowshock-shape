"""
The final stage: correct the *-forma.reg ds9
region files
Please run at the PC-orig folder
"""


import os
import glob
import json
import astropy.coordinates as coords
import astropy.units as u

data = json.load(open("PC_offsets.json"))
path_to_correct = "../PC-correct"


def has_region(line):
    """
    True if line looks like a ds9 region
    """
    return "(" in line and ")" in line

def correct_coords(line,d):
    """
    Apply correction in a ds9 file
    """
    ra,dec = line.split(",")[:2]
    ra = ra.split("(")[-1]
    dec = dec.split(")")[0]
    c = coords.ICRSCoordinates(ra,dec,unit=(u.hour,u.deg))
    new_ra = (c.ra-d["avdx"][0]*u.arcsec).to_string(sep=":")
    new_dec = (c.dec-d["avdy"][0]*u.arcsec).to_string(sep=":")
    return line.replace(ra,new_ra).replace(dec,new_dec)

for field,d in data.items():
    formas = glob.glob(os.path.join(field, "*-forma.reg"))
    masks = glob.glob(os.path.join(field, "*-mask.reg"))
    regfiles = formas + masks
    wcsfield = os.path.join(path_to_correct,field)
    if not os.path.isdir(wcsfield):
        print "*****Creating {} directory".format(wcsfield)
        os.mkdir(wcsfield)
    for regfile in regfiles:
        newfile = os.path.join(wcsfield,regfile.split("/")[-1])
        with open(regfile) as f:
            lines = f.readlines()
        newlines = []
        for line in lines:
            # Corrected regions have new colors
            if "green" in line:
                line = line.replace("green","yellow")
            if "magenta" in line:
                line = line.replace("magenta","red")
            if has_region(line):  
                newlines.append(correct_coords(line,d))
            else:
                newlines.append(line)
        print "****writting {}****".format(newfile)
        with open(newfile,"w") as f:
            f.writelines(newlines)



