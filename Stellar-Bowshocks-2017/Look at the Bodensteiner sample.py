# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import numpy as np
import pandas as pd
from scipy import stats
from matplotlib import pyplot as plt
import seaborn as sns
from astropy.table import Table, join
from pathlib import Path
from astropy.io import fits

sns.set_color_codes()
sns.set_context("poster")
pd.set_option('display.precision', 3)
pd.set_option("display.max_columns", 999)
pd.set_option("display.max_rows", 9999)
# -

# Get the tables directly from Vizier is best.  This has the advantage that the RA and Dec are in a single column each.  Although the disadvantage is that all the string columns are now bytes, not unicode.

from astroquery.vizier import Vizier

Vizier.ROW_LIMIT = -1
catalogs = Vizier.get_catalogs("J/A+A/618/A110")

catalogs

# There are two tables, which cover the same 255 sources but with different data.

catalogs[1]

catalogs[1].dtype

catalogs[0]

# Make a mask to select the bow shock (bs) and non-aligned bow shock (bsna) sources.  There are 94 of them.

table = join(catalogs[0], catalogs[1])

table.sort(keys=["RAJ2000", "DEJ2000"])

m = (table["MClass"] == "bs") | (table["MClass"] == "bsna")
m.sum()

# Now make a table that merges the two tables but only for the selected 94 sources.

table = table[m]

table

import astropy.units as u
import astropy.coordinates as coord


def skycoord_from_table_row(data):
    ra = data["RAJ2000"]
    dec = data["DEJ2000"]
    return coord.SkyCoord(f'{ra} {dec}', unit=(u.hourangle, u.deg))



skycoord_from_table_row(table[-1])

# Try out obtaining the WISE images from IPAC

# +
import requests
import io

SIA_URL = 'https://irsa.ipac.caltech.edu/SIA'

sia_params = {
    'COLLECTION': 'wise_allwise',
    'RESPONSEFORMAT': 'VOTABLE',
    'FORMAT': 'image/fits',
    'POS': 'circle $RA $DEC 0.0',
}
# -

source_data = table[-1]
c = skycoord_from_table_row(source_data)
sia_params['POS'] = f"circle {c.to_string()} 0.0"


# +
r = requests.get(SIA_URL, params=sia_params)

tab = Table.read(io.BytesIO(r.content), format='votable')

# -

tab

# Expand the image size for bigger bows
expand = 1.0
for threshold in 40.0, 80.0, 160.0:
    if 60*source_data["R0A"] > threshold:
        expand *= 2

expand

source_data["R0A"]

BASE_IMAGE_SIZE_ARCMIN = 8.0

image_size = BASE_IMAGE_SIZE_ARCMIN*expand
image_params = {
    "center": f"{c.ra.deg:.4f},{c.dec.deg:.4f}",
    "size": f"{image_size}, {image_size} arcmin",
    "gzip": 0,
}

OUTPUT_IMAGE_DIR = Path(".")

f"{source_data['HD']:06d}"

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

table[table["HD"] == 111123]

table[table["HD"] == 115842]

table[table["HD"] == 213087]

table[table["HD"] == 214680]

table[table["HD"] == 203064]

table[table["HD"] == 175362]

table[table["HD"] == 34078]


