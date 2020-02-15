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

from astropy.table import Table, join
from pathlib import Path
import numpy as np

SOURCE_DIR = Path('OB/Kobulnicky2016')
source_table = Table.read(
    str(SOURCE_DIR / 'table1.dat'),
    format='ascii.cds',
    readme=str(SOURCE_DIR / 'ReadMe')
)

star_table = Table.read('star-ratings.tab', format='ascii.tab')
combo_table = join(source_table, star_table, join_type='outer')

# Look for big ones, but only those that have been rated.  The others are WISE-only sources, which we aren't using yet

mbig = combo_table['R0'] > 40.0
mrated = ~combo_table['Rating'].mask

# I have put the cut-off at 40 instead of 60 because even some of the 40 arcsec ones would benefit from a bit more room to trace the wings.

combo_table[mbig & mrated]["Seq", "Name", "Alias", "R0", "Rating", "Comment"]

# Many of these are missing. This was because they are only in WISE.  Now that I have all the WISE images too, we can look at the new big ones:

combo_table[mbig & ~mrated]["Seq", "Name", "Alias", "R0", "Rating", "Comment"]


