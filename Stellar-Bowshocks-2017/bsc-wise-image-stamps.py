import sys
import io
from pathlib import Path
import requests
import numpy as np
from astropy.table import Table, join
from astropy.io import fits
import astropy.units as u
import astropy.coordinates as coord
from astroquery.vizier import Vizier

SIA_URL = 'https://irsa.ipac.caltech.edu/SIA'

sia_params = {
    'COLLECTION': 'wise_allwise',
    'RESPONSEFORMAT': 'VOTABLE',
    'FORMAT': 'image/fits',
    'POS': 'circle $RA $DEC 0.0',
}

Vizier.ROW_LIMIT = -1
catalogs = Vizier.get_catalogs("J/A+A/618/A110")
source_table = join(catalogs[0], catalogs[1])
source_table.sort(keys=["RAJ2000", "DEJ2000"])
# Restrict to only bow shock sources
m = (source_table["MClass"] == "bs") | (source_table["MClass"] == "bsna")
source_table = source_table[m]

OUTPUT_IMAGE_DIR = Path('OB/BSC-WISE')
OUTPUT_IMAGE_DIR.mkdir(exist_ok=True)

BASE_IMAGE_SIZE_ARCMIN = 8.0

def skycoord_from_table_row(data):
    ra = data["RAJ2000"]
    dec = data["DEJ2000"]
    return coord.SkyCoord(f'{ra} {dec}', unit=(u.hourangle, u.deg))

try:
    k1 = int(sys.argv[1])
except:
    k1 = 1
try:
    k2 = int(sys.argv[2])
except:
    k2 = None

# Loop over all sources in the table
for source_data in source_table[k1-1:k2]:
    print(source_data["HD", "Name", "R0A"])

    # Make a SkyCoord object
    c = skycoord_from_table_row(source_data)
    sia_params['POS'] = f"circle {c.to_string()} 0.0"

    # Perform a search around the specified coordinates
    r = requests.get(SIA_URL, params=sia_params)

    tab = Table.read(io.BytesIO(r.content), format='votable')

    # Expand the image size for bigger bows
    expand = 1.0
    for threshold in 40.0, 80.0, 160.0:
        if 60*source_data["R0A"] > threshold:
            expand *= 2
    image_size = BASE_IMAGE_SIZE_ARCMIN*expand
    image_params = {
        "center": f"{c.ra.deg:.4f},{c.dec.deg:.4f}",
        "size": f"{image_size}, {image_size} arcmin",
        "gzip": 0,
    }
    # Now fetch images in each band
    for data in tab:
        print(
            f"Fetching image ({image_size} arcmin square) from",
            data['access_url'].decode(),
        )
        r = requests.get(data['access_url'], params=image_params)
        hdulist = fits.open(io.BytesIO(r.content))
        # Get name of WISE bandpass as a unicode string
        bpname = data['energy_bandpassname'].decode()
        hdulist.writeto(
            OUTPUT_IMAGE_DIR / f"HD{source_data['HD']:06d}-{bpname}.fits",
            overwrite=True,
        )
