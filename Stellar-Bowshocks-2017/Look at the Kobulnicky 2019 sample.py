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
from astroquery.vizier import Vizier

sns.set_color_codes()
sns.set_context("poster")
pd.set_option('display.precision', 3)
pd.set_option("display.max_columns", 999)
pd.set_option("display.max_rows", 9999)
# -

Vizier.ROW_LIMIT = -1
catalogs = Vizier.get_catalogs("J/AJ/158/73")
catalogs

catalogs[0]


