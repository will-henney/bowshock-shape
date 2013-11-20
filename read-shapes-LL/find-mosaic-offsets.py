"""
The main goal is to find the position offsets
in the GO*** using the Robberto coordinates as
reference, using the same method as the 
find-bally-offsets.py program
"""
import os
import numpy as np
import glob
import json
from astropy import units as u
from astropy import coordinates as coord
import image_statistics


def coords(data):
    """
    Extract RA and Dec from data for a single source
    """
    ra = data["RA"]
    dec = data["Dec"]
    return coord.ICRSCoordinates(ra,dec,unit = (u.hour,u.degree))


# All the mosaics have the same offset, so, it is only 
# neccesary to find the offset in one image


# Find the coordinates in the "shifted" system using the proplyds
# JSON files

json_patterns = "*-arcdata.json"
table_file = "ll-data.json"
#1 Be sure of read correctly the table and find all the JSON files

# table_file path   path_to/Dropbox/ProplydMIR/LuisLL
# The program must be run in the folder path_to/Dropbox/ProplydMIR/Jorge_prop/
table_path = os.path.join("../LuisLL",table_file)
table = json.load(open(table_path))

# print glob.glob(json_patterns) this command works wonderfully =)
dx = []
dy = []
sources = []


for sourcefile in glob.glob(json_patterns):
    source = sourcefile.replace("-arcdata.json","")
    #print source so far working
    source_db = json.load(open(sourcefile))
    pcoords = coords(source_db["star"])
    print source,pcoords

    if source in table:
        rcoords = coords(table[source])
        dRa = pcoords.ra - rcoords.ra
        dDec = pcoords.dec - rcoords.dec
        print "Offset:", dRa,dDec
        dx.append(dRa)
        dy.append(dDec)
        sources.append(source)
#print dx,dy
avdx, sigdx = image_statistics.trimean_and_iqr(dx)
avdy, sigdy = image_statistics.trimean_and_iqr(dy)

#Apparently dx and dy cannot be converted into arrays, so i'm having troubles

#offsets = {
#    "dx":dx,
#    "dy":dy,
#    "sources":sources,
#    "avdx": [avdx,sigdx],
#    "avdy": [avdy,sigdy]
#    }
#with open("mosaic_offsets.json",w) as f:
#    json.dump(offsets,f,indent=2)

