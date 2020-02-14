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

import os
import sys
import io
import requests
import xmltodict
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import astropy.coordinates as coord

# +
SIA_URL = 'https://irsa.ipac.caltech.edu/SIA'

sia_params = {
    'COLLECTION': 'wise_allwise',
    'RESPONSEFORMAT': 'VOTABLE',
    'FORMAT': 'image/fits'
}
# -

SOURCE_DIR = 'OB/Kobulnicky2016'
source_table = Table.read(
    os.path.join(SOURCE_DIR, 'table1.dat'),
    format='ascii.cds',
    readme=os.path.join(SOURCE_DIR, 'ReadMe')
)

OUTPUT_IMAGE_DIR = 'OB/WISE'
IMAGE_SIZE_DEGREES = 4.0/60.0           


def skycoord_from_table_row(data):
    ra = f"{data['RAh']} {data['RAm']} {data['RAs']}"
    dec = f"{data['DE-']}{data['DEd']} {data['DEm']} {data['DEs']}"
    return coord.SkyCoord(f'{ra} {dec}', unit=(u.hourangle, u.deg))


# Use zeta Oph as an example. 

source_data = source_table[12]

source_data

# Make a SkyCoord object for specifying location
c = skycoord_from_table_row(source_data)
f"circle {c.to_string()} 0.0"

# Perform a search around the specified coordinates
r = requests.get(SIA_URL, 
                 params={**sia_params, 'POS': f"circle {c.to_string()} 0.0"})

tab = Table.read(io.BytesIO(r.content), format='votable')

tab


