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

# The purpose of this file is twofold:
#
# 1. Calculate derived quantities such as planitude, alatude, kurtosis, etc
# 2. Combine the results from the different circle-fit variants ("peak" vs "both" for finding the ridge, and 50, 60, 60 for `THETA_CIRC`)
#
# From the latter, we will also calculate better estimates of parameter uncertainty that include the variation over the fit variants. 

import numpy as np
import pandas as pd
from scipy import stats

from astropy.table import Table

# These are the different variants of the circle fit: the maximum $\theta$ from the axis that is included in the fit, and whether the ridge is identified using just the peak radius, or the peak and the mean radius. 

circle_variants = ["50", "60", "70"]
ridge_variants = ["peak", "both"]

# Read the data files as `astropy.table`s and convert to pandas. 

data = {}
for circle in circle_variants:
    for ridge in ridge_variants:
        tab = Table.read(
            f"mipsgal-arcfit-0{circle}_{ridge}.tab",
            format='ascii.tab',
        )
        data[(circle, ridge)] = tab.to_pandas()

# Make a long-form dataset of the whole lot by `concat`ing the individual tables, which makes a multi-index, and then doing `reset_index` that breaks out the levels of the multi-index as new columns, which we call `DTHC` and `RIDGE`.  We also take the opportunity to drop the coordinate columns and other Kobulnicky data that we are not using here. 

pd.set_option('display.precision', 3)
pd.set_option("display.max_columns", 999)

unwanted = ["RAh", "RAm", "RAs", "DE-", "DEd", "DEm", "DEs", ]
unwanted += ["Ref", "Alias", "8um", "Unc", "R0", "Hmag", "4.5mag", "Ak", "Env", "Name", "PA"]
df = pd.concat(data, names=['DTHC', 'RIDGE']).reset_index().drop(
    columns=["level_2"] + unwanted,
)
df


# Helper function to do log transform of a quantity and its uncertainty, and also for the ratio of two columns and their uncertainties.

# +
def logify(x, dx):
    """Base-10 log of quantity `x` and its error `dx`"""
    return np.log10(x), np.log10(np.e)*dx/x

def logratio(x, dx, y, dy):
    """Base-10 log of ratio (`x` +/-`dx`) / (`y` +/-`dy`) plus error"""
    r = x/y
    dr_over_r = np.sqrt((dx/x)**2 + (dy/y)**2)
    return np.log10(r), np.log10(np.e)*dr_over_r

logratio(2.0, 0.5, 2.0, 0.5)
# -

# Test that the function works for adding new columns to a dataframe.

d = pd.DataFrame.from_dict(dict(
    A=[1, 2, 3], dA=[0.1, 0.2, 0.3],
    B=[3, 2, 1], dB=[0.1, 0.2, 0.3]
))
d["log R"], d["d log R"] = logratio(d["A"], d["dA"], d["B"], d["dB"])
d

# For the various varieties of errors in a quantity `X`, we will use the following naming convention:
#
# * `e X` is the error from the `_sigma` entries, which mainly come from the variations in the ridge radius over $\pm 10^\circ$ around 0 or 90 degrees.
# * `d X` is the asymmetry between the two sides (e.g., for R90)
# * `E X` is the dispersion over the different circle-fit variants
#
# We calculate these for planitude and alatude. 

df["log Pi"], df["e log Pi"] = logratio(
    df["Rc"], df["Rc_sigma"],
    df["R0_fit"], df["R0_sigma"]
) 
df["log Lam"], df["e log Lam"] = logratio(
    0.5*(df["R90p"] + df["R90n"]), 0.5*(df["R90p_sigma"] + df["R90n_sigma"]),
    df["R0_fit"], df["R0_sigma"]
) 
df["d log Lam"] = np.log10(1 + 0.5*(df["R90p"] - df["R90n"])/df["R0_fit"])
df["log R0"], df["d log R0"] = logify(df["R0_fit"], df["R0_sigma"])

# And for the shell thickness.  We use the minimum of the direct and Gaussian methods, with uncertainty given by the absolute difference.

H = 0.5*(df["H_g"] + df["H_d"])
eH = 0.5*np.abs(df["H_g"] - df["H_d"])
df["log H"], df["e log H"] = logratio(
    H, eH,
    df["R0_fit"], df["R0_sigma"]
) 

# For the angular breadth, we don't take log for now (**should we?**).  We also calculate ratios of the widths at different levels, which indicate kurtosis:
#
# * Platykurtic (flat-topped): core > 0.5 and wing < 1.5
# * Leptokurtic (peak and wings): core < 0.5 and wing > 1.5

# +
df["breadth50"] = 0.5*(df["th_p"] - df["th_m"])
df["d_breadth50"] = 0.5*(df["th_p"] + df["th_m"])

df["core_bratio"] = 0.5*(df["th_p_841"] - df["th_m_841"])/df["breadth50"]
df["wing_bratio"] = 0.5*(df["th_p_210"] - df["th_m_210"])/df["breadth50"]
# -

# Finally, we eliminate the columns that we no longer need.

df = df.drop(
    columns=[
        "Rc", "Rc_sigma", "R90p", "R90n", "R90p_sigma", "R90n_sigma",
        "R0_fit", "R0_sigma", "pa_circ", "delta_pa", 
        "th_p", "th_m", "th_p_841", "th_m_841", "th_p_210", "th_m_210", 
        "th_g", "dth_g", "H_g", "H_d", "R0_g",
    ]
)
df

# We can use `groupby` on the fit variant columns to see the systematic trends with these variants.  We use `query` to restrict the results to 3, 4, 5-star sources, and we use `dropna` to eliminate rows with NaNs, but only checking the most important columns:

variant_groups = (
    df.query('Rating >=3')
    .dropna(subset=["log Pi", "log Lam", "log R0"])
    .drop(columns=['Seq', 'Rating'])
    .groupby(['RIDGE', 'DTHC'])
)
variant_groups.agg([np.mean, np.std, 'count'])

variant_groups.agg([np.mean, np.std]).plot(
    y=[("log Pi", "mean"), ("log Lam", "mean")],
    use_index=True,
    kind="bar",
    ylim=[0.10, 0.40],
)

variant_groups.agg([np.median, stats.iqr])

variant_groups.agg([np.median, stats.iqr]).plot(
    y=[("log Pi", "median"), ("log Lam", "median")],
    use_index=True,
    kind="bar",
    ylim=[0.10, 0.40],
)

# We did two versions, one with mean and stddev, the other with median and iinterquartile range.  Results are shown in tables and are plotted for planitude and alatude.  Strangely, the ridge strategy makes more difference to planitude than the $\theta$ range, especially for the mean (for the median, they both have an effect).

pd.set_option("display.max_rows", 999)


# The next thing we can do is `.groupby(['Seq'])`, so we can get the average and dispersion over the fit variants for each source:

# +
def compress_double_levels(colnames):
    return colnames[0] if "mean" in colnames[1] else "E " + colnames[0]
    
ddf = df.query("Rating >= 3").groupby(['Seq']).agg([np.mean, np.std])
ddf.columns = ddf.columns.map(compress_double_levels)
ddf
# -

# We then add another `groupby`, by `Rating`.  This time, we allow all ratings, but we drop NaNs for some columns. 

ddf = (
    df
    .dropna(subset=["log Pi", "log Lam", "log R0"])
    .groupby(['Seq'])
    .agg([np.mean, np.std])
    .drop(columns=[('Rating', 'std'), ])
)
ddf.columns = ddf.columns.map(compress_double_levels)
d_by_star = (
    ddf
    .groupby(['Rating'])
    .agg([np.mean, np.std])
)
d_by_star

# From this table, we can see the following:
#
# * Star-rating is positively correlated with `Peak24`, `Contrast`, `log Pi`, `log Lam`, `log R0`, `core_bratio`
# * Negative corelations are mainly with the `E` and `e` dispersions, but are also with the population stddev of `log Lam`, `log Pi`, `core bratio`, and `breadth50`
#
#

pd.DataFrame(d_by_star[1:].reset_index().corr()['Rating']).query("Rating > 0.8")

pd.DataFrame(d_by_star.reset_index().corr()['Rating']).query("Rating < -0.8")

# Go back to 3, 4, 5 star and look at correlations.   Also, do some sigma clipping on the errors.

import seaborn as sns
sns.set_color_codes()
sns.set_context("poster")

ddf = (
    df
    .query('Rating >= 0')
    .dropna(subset=["log Pi", "log Lam", "log R0"])
    .groupby(['Seq'])
    .agg([np.mean, np.std])
    .drop(columns=[('Rating', 'std'), ])
)
ddf.columns = ddf.columns.map(compress_double_levels)
ddf = ddf.query(
    '`d log R0` < 0.1 & `E log Pi` < 0.35 & `E log Lam` < 0.1 & `E log H` < 0.2 & `e log H` < 0.2'
)
ddf.describe()

cols = [
    "log Pi", "log Lam", "log R0",
    "E log Pi", "E log Lam", "E log R0",
    "e log Pi", "e log Lam", "d log R0",
    "e log H", "E log H",
]
sns.pairplot(
    ddf, 
    hue="Rating", 
    vars=cols, 
    palette="magma",
    plot_kws=dict(s=50), 
    diag_kind="hist",
)

# This makes it look like it might be possible to put a cut on `E log Pi` < 0.2, perhaps.  Or even < 0.1. `

qstrict = " and ".join([
    "`E log Pi` < 0.1",
    "`E log Lam` < 0.05",
    "`E log R0` < 0.05",
    "`e log Pi` < 0.05",
    "`e log Lam` < 0.1",
    "`e log H` < 0.15",
    "`E log H` < 0.05",
    "`log H` > -1 ",
    "wing_bratio < 2",
])

ddf.query(qstrict).describe()

sns.pairplot(
    ddf.query(qstrict), 
    hue="Rating", 
    vars=cols, 
    palette="magma",
    plot_kws=dict(s=50), 
    diag_kind="hist",
)

sns.pairplot(
    ddf.query("`E log Pi` < 0.1 & `log H` > -1"), 
    hue="Rating", 
    vars=["log Pi", "log Lam", "log H", "breadth50"], 
    plot_kws=dict(s=30),
    palette="magma",
    diag_kind="hist",
)

# Look at the sources that have small shell widths.  Some of these are actually OK, such as K292.

df.query("`log H` < -0.5").sort_values(['Seq'])

# This is an example of a source that goes badly wrong on the peak fits.  As a result, it fails the `E log H` test and so is excluded from `ddf`.

df.query('Seq == 587')

ddf[["log Pi", "log Lam", "log H", "breadth50"]].corr()

ddf.query("`E log Pi` < 0.1 & `log H` > -1")[["log Pi", "log Lam", "log H", "breadth50"]].corr()

ddf.query(qstrict)[["log Pi", "log Lam", "log H", "breadth50"]].corr()

sns.pairplot(
    ddf.query("`E log Pi` < 0.1 & `log H` > -1 & wing_bratio < 2"), 
    hue="Rating", 
    vars=["breadth50", "core_bratio", "wing_bratio", "log Lam"], 
    diag_kind="hist",
    palette="magma",
    plot_kws=dict(s=30)
)

ddf.query(
    "`E log Pi` < 0.1 & `log H` > -1 & wing_bratio < 2"
)[["log Pi", "log Lam", "breadth50", "core_bratio", "wing_bratio"]].corr()

ddf.query(
    qstrict
)[["log Pi", "log Lam", "breadth50", "core_bratio", "wing_bratio"]].corr()

sns.pairplot(
    ddf.query(qstrict), 
    hue="Rating", 
    vars=["log Pi", "log Lam", "log H", "breadth50"], 
    palette="magma",
    plot_kws=dict(s=30),
    diag_kind="hist",
)

sns.pairplot(
    ddf.query(qstrict), 
    hue="Rating", 
    vars=["breadth50", "core_bratio", "wing_bratio", "log Lam"], 
    diag_kind="hist",
    palette="magma",
    plot_kws=dict(s=30)
)

# So the moral of this is that it might not be worth imposing the strict cut-offs.  On the other hand, they are a viable alternative to excluding 1, 2-star sources. 

# ## Write the result to a file
#
# We should save it as an astropy Table for compatibility with before.

Table.from_pandas(ddf).write("mipsgal-arcfit-all-variants.tab", format='ascii.tab', overwrite=True)


