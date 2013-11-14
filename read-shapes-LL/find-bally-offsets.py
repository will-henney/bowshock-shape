"""Find the offsets between the Bally images and the 'true' positions

We do this by reading in all the source '*-xy.json' files that we
can find and extracting the coordinates of the source. 

"""
import glob
import os
import json
from astropy import units as u 
from astropy import coordinates as coord
import image_statistics

def coords(data):
    "Extract star RA and Dec from database for a single source"
    ra = data["RA"]
    dec = data["Dec"]
    return coord.ICRSCoordinates(ra, dec, unit=(u.hour, u.degree))


BALLY_FILE_PATTERN = "j8oc??010_drz"
DB_FILE_PATTERN = "*-xy.json"
TABLE_FILE = "ll-data.json"

if __name__ == "__main__":

    # First read in the table
    table = json.load(open(TABLE_FILE))

    offsets = {}
    for field in glob.glob(BALLY_FILE_PATTERN):
        bally_id = field[4:6]
        offsets[bally_id] = {"file": field}
        print "**************** ", field, " ****************" 

        dx = []; dy = []; sources = []
        for sourcefile in glob.glob(
                os.path.join(field, DB_FILE_PATTERN)
        ):
            source = sourcefile.split('/')[-1].replace("-xy.json", "")
            source_db = json.load(open(sourcefile))
            bcoords = coords(source_db["star"])
            print source, bcoords
            if source in table:
                rcoords = coords(table[source])
                dRA = bcoords.ra - rcoords.ra
                dDec = bcoords.dec - rcoords.dec
                print "Offset: ", dRA.arcsec, dDec.arcsec
                dx.append(dRA.arcsec)
                dy.append(dDec.arcsec)
                sources.append(source)
            else:
                print "NOT FOUND IN TABLE!"

        avdx, sigdx = image_statistics.trimean_and_iqr(dx)
        avdy, sigdy = image_statistics.trimean_and_iqr(dy)

        offsets[bally_id].update(
            dx=dx, dy=dy, sources=sources,
            avdx=[avdx, sigdx], avdy=[avdy, sigdy]
        )
            

    with open("bally-offsets.json", "w") as f:
        json.dump(offsets, f, indent=2)
