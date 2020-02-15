import sys
import io
from pathlib import Path
import requests
import numpy as np
from astropy.table import Table
from astropy.io import fits
import astropy.units as u
import astropy.coordinates as coord

SIA_URL = 'https://irsa.ipac.caltech.edu/SIA'

sia_params = {
    'COLLECTION': 'wise_allwise',
    'RESPONSEFORMAT': 'VOTABLE',
    'FORMAT': 'image/fits',
    'POS': 'circle $RA $DEC 0.0',
}

SOURCE_DIR = Path('OB/Kobulnicky2016')
source_table = Table.read(
    str(SOURCE_DIR / 'table1.dat'),
    format='ascii.cds',
    readme=str(SOURCE_DIR / 'ReadMe')
)

OUTPUT_IMAGE_DIR = Path('OB/WISE')
OUTPUT_IMAGE_DIR.mkdir(exist_ok=True)

BASE_IMAGE_SIZE_ARCMIN = 8.0

def skycoord_from_table_row(data):
    ra = f"{data['RAh']} {data['RAm']} {data['RAs']}"
    dec = f"{data['DE-']}{data['DEd']} {data['DEm']} {data['DEs']}"
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
    print(source_data["Seq", "Name", "Alias", "R0"])

    # Make a SkyCoord object
    c = skycoord_from_table_row(source_data)
    sia_params['POS'] = f"circle {c.to_string()} 0.0"

    # Perform a search around the specified coordinates
    r = requests.get(SIA_URL, params=sia_params)

    tab = Table.read(io.BytesIO(r.content), format='votable')

    # Expand the image size for bigger bows
    expand = 1.0
    for threshold in 40.0, 80.0, 160.0:
        if source_data["R0"] > threshold:
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
            OUTPUT_IMAGE_DIR / f"{source_data['Seq']:04d}-{bpname}.fits",
            overwrite=True,
        )
