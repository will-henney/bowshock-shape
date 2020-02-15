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
#
# *Update: 2019-02-19* We no longer use the `both` variant, just the `peak`.

circle_variants = ["50", "60", "70"]

# Read the data files as `astropy.table`s and convert to pandas. 

data = {}
for circle in circle_variants:
    tab = Table.read(
        f"mipsgal-arcfit-0{circle}_peak.tab",
        format='ascii.tab',
    )
    data[circle] = tab.to_pandas()

# Make a long-form dataset of the whole lot by `concat`ing the individual tables, which makes a multi-index, and then doing `reset_index` that breaks out the levels of the multi-index as new columns, which we call `DTHC` ~and `RIDGE`~.  We also take the opportunity to drop the coordinate columns and other Kobulnicky data that we are not using here. 

pd.set_option('display.precision', 3)
pd.set_option("display.max_columns", 999)

unwanted = ["RAh", "RAm", "RAs", "DE-", "DEd", "DEm", "DEs", ]
unwanted += ["Ref", "Alias", "8um", "Unc", "R0", "Hmag", "4.5mag", "Ak", "Env", "Name", "PA"]
df = pd.concat(data, names=['DTHC']).reset_index().drop(
    columns=["level_1"] + unwanted,
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

# And for the shell thickness.  We use outer half-thickness only, since the inner one is very sensitive to instrumental resolution.

# Temporary fix to the half-widths - REMOVE once everything is re-run
df["H_d_out"] += df["R0_fit"]
df["H_d_in"] -= df["R0_fit"]

# Check the half-thickness correlation

df.query("Rating >= 4")[["H_d_out", "H_d_in"]].corr()

df.query("Rating >= 4")[["H_d_out", "H_d_in"]].describe()

# It is a bit annoying that there is so little differenc between these.

H = df["H_d_out"]
eH = 0.5*np.abs(0.5*df["H_g"] - H)
df["log H"], df["e log H"] = logratio(
    H, eH,
    df["R0_fit"], df["R0_sigma"]
) 

# Also add linear versions for later
df["log h"], df["e log h"] = logify(H, eH)

# And, we also save a version of the inner half-width, just in case.

Hin = np.abs(df["H_d_in"])
eHin = 0.5*np.abs(0.5*df["H_g"] - Hin)
df["log Hin"], df["e log Hin"] = logratio(
    Hin, eHin,
    df["R0_fit"], df["R0_sigma"]
) 
df["log hin"], df["e log hin"] = logify(Hin, eHin)


# For the angular breadth, we don't take log for now (**should we?**).  We also calculate ratios of the widths at different levels, which indicate kurtosis:
#
# * Platykurtic (flat-topped): core > 0.5 and wing < 1.5
# * Leptokurtic (peak and wings): core < 0.5 and wing > 1.5

# **New:** Define a function to convert from angular breadth to arc length.  Arclength $s = R_c \theta_c$, where $R_c \sin\theta_c = R(\theta) \sin\theta$ and 
# $R_c - R_0 + R(\theta) \cos\theta = R_c \cos\theta_c$.
#
# This has the solution:
# $$
# s = \Pi \arccos\left[ 
# \cos\theta \left(1 - p^2 \sin^2\theta\right)^{1/2}
# + p \sin^2\theta
# \right]
# $$
# where $p = 1 - \Pi^{-1}$.
#

def arclength_from_angle(d, theta, Rc, R0):
    t = np.deg2rad(d[theta])
    ct, st = np.cos(t), np.sin(t)
    Pi = d[Rc]/d[R0]
    p = 1.0 - 1.0/Pi
    muc = ct*np.sqrt(1.0 - p**2*st**2) + p*st**2
    return  np.sign(t)*d[Rc]*np.arccos(muc)



# And there is an alternative version that uses the alatude instead of planitude, and which assumes a parabola shape. 

def arclength_from_angle_parabola(d, theta, R90, R0):
    th = np.deg2rad(d[theta])
    cott = np.cos(th)/np.sin(np.abs(th))
    Lam = d[R90]/d[R0]
    t = (2./Lam)*np.sqrt(1 + (Lam*cott/2)**2) - cott
    t *= np.sign(th)
    s = (Lam/2)**2 * (t*np.sqrt(1 + t**2) + np.arcsinh(t))
    return d[R0]*s



d = pd.DataFrame({
    "R0": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    "Rc": [1.0, 1.0, 3.0, 1.0, 1.0, 1.8, 5.0, 5.0],
    "R90": [1.0, 1.0, 2.0, 1.5, 1.5, 1.3, 3.0, 3.0],
    "th": [30, -30, 90, 170.0, -170.0, 134.0, 90.0, 170.0],
})
d["s"] = arclength_from_angle(d, "th", "Rc", "R0")
d["sp"] = arclength_from_angle_parabola(d, "th", "R90", "R0")
d["ss"] = d["R0"]*np.deg2rad(d["th"])
d

# This is consistent with the table that I calculated in `stellar-bowshocks.org`

# Use it to find half arc-lengths in physical units

# +
s_p = arclength_from_angle(df, "th_p", "Rc", "R0_fit")
s_m = arclength_from_angle(df, "th_m", "Rc", "R0_fit")
s_p_841 = arclength_from_angle(df, "th_p_841", "Rc", "R0_fit")
s_m_841 = arclength_from_angle(df, "th_m_841", "Rc", "R0_fit")
s_p_210 = arclength_from_angle(df, "th_p_210", "Rc", "R0_fit")
s_m_210 = arclength_from_angle(df, "th_m_210", "Rc", "R0_fit")

ss_p = arclength_from_angle_parabola(df, "th_p", "R90p", "R0_fit")
ss_m = arclength_from_angle_parabola(df, "th_m", "R90n", "R0_fit")
ss_p_841 = arclength_from_angle_parabola(df, "th_p_841", "R90p", "R0_fit")
ss_m_841 = arclength_from_angle_parabola(df, "th_m_841", "R90n", "R0_fit")
ss_p_210 = arclength_from_angle_parabola(df, "th_p_210", "R90p", "R0_fit")
ss_m_210 = arclength_from_angle_parabola(df, "th_m_210", "R90n", "R0_fit")
# -

# Save the angular breadth as before.

df["breadth50"] = 0.5*(df["th_p"] - df["th_m"])
df["d_breadth50"] = 0.5*(df["th_p"] + df["th_m"])

# Average half arc-length and estimate of its error, plus log versions of both the physical breadth and its ratio to $R_0$.  Note that the `s` version is using the circle/planitude arc-length estimate, while the `ss` version is using the parabola/alatude arc-length estimate.

# +
s = 0.5*(s_p - s_m)
es = 0.5*(s_p + s_m)
df["log s"], df["e log s"] = logify(s, es)
df["log S"], df["e log S"] = logratio(
    s, es,
    df["R0_fit"], df["R0_sigma"]
) 

ss = 0.5*(ss_p - ss_m)
ess = 0.5*(ss_p + ss_m)
df["log ss"], df["e log ss"] = logify(ss, ess)
df["log SS"], df["e log SS"] = logratio(
    ss, ess,
    df["R0_fit"], df["R0_sigma"]
) 

# -

df.query("Rating >= 3")[["log S", "log SS"]].corr()

# It is interesting that the two breadth measurements are *much better correlated than for WISE*

df.query("Rating >= 3")[["log S", "log SS"]].describe()

# The kurtosis is now defined in terms of the physical lengths.  This needs revamping to make it log and to use both arc-length estimates, but that is too much work for now.

df["core_bratio"] = 0.5*(s_p_841 - s_m_841)/s
df["wing_bratio"] = 0.5*(s_p_210 - s_m_210)/s

# Finally, we eliminate the columns that we no longer need.

df = df.drop(
    columns=[
        "Rc", "Rc_sigma", "R90p", "R90n", "R90p_sigma", "R90n_sigma",
        "R0_fit", "R0_sigma", "pa_circ", "delta_pa", 
        "th_p", "th_m", "th_p_841", "th_m_841", "th_p_210", "th_m_210", 
        "th_g", "dth_g", "H_g", "H_d", "R0_g", "H_d_out", "H_d_in",
    ]
)

df

# We can use `groupby` on the fit variant columns to see the systematic trends with these variants.  We use `query` to restrict the results to 3, 4, 5-star sources, and we use `dropna` to eliminate rows with NaNs, but only checking the most important columns:

variant_groups = (
    df.query('Rating >=3')
    .dropna(subset=["log Pi", "log Lam", "log R0"])
    .drop(columns=['Seq', 'Rating'])
    .groupby('DTHC')
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
# * Star-rating is highly positively correlated with `Peak24`, ~`log Pi`~, `log Lam`, `log R0`, ~`core_bratio`~ plus (new for WISE) `log h`
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
ddf.query('Rating >= 3').describe()

cols = [
    "log Pi", "log Lam", "log R0",
    "E log Pi", 
    "e log Pi", "e log Lam", "d log R0",
    "e log H", 
]
sns.pairplot(
    ddf, 
    hue="Rating", 
    vars=cols, 
    palette="magma",
    plot_kws=dict(s=50), 
    diag_kind="hist",
)

# This makes it look like it might be possible to put a cut on `E log Pi` < 0.2, perhaps.  Or even < 0.1.

ddf[["log Pi", "log Lam", "log H", "log S", "log SS"]].corr()

ddf.query(
    "`E log Pi` < 0.1 & `log H` > -1 & `log S` > -0.3 & Rating >= 3"
)[["log Pi", "log Lam", "log H", "log S", "log SS"]].corr()

sns.pairplot(
    ddf.query("Rating >= 3").query(
        "`E log Pi` < 0.1 & `log H` > -1 & wing_bratio < 2 & -0.1 < `log S` < 0.8 & -0.1 < `log SS` < 0.8"
    ), 
    hue="Rating", 
    vars=["log S", "log SS", "breadth50", "core_bratio", "wing_bratio", "log Pi", "log H"], 
    diag_kind="kde",
    palette="magma_r",
    plot_kws=dict(s=30)
)

# ## Write the result to a file
#
# We should save it as an astropy Table for compatibility with before.

Table.from_pandas(ddf.reset_index()).write(
    "mipsgal-arcfit-peak-variants.tab", 
    format='ascii.tab', 
    overwrite=True
)

# Note `astropy.table.Table.from_pandas` ignores the index, so that we had to do `reset_index` to get `Seq` o be an actual column.

ddf.reset_index().head()

# ## What if we were to only use the peak method?
#
# This will reduce the `E` values, and will shift $\Pi$ and $\Lambda$ to larger values.
#
#

# But this is now done from the start, however we still re-create the dataframe here, since we want to undo the cut-offs on erros that we had applied above.

dfp = (
    df
    .query('Rating >= 0')
    .dropna(subset=["log Pi", "log Lam", "log R0"])
    .groupby(['Seq'])
    .agg([np.mean, np.std])
    .drop(columns=[('Rating', 'std'), ])
)
dfp.columns = dfp.columns.map(compress_double_levels)

select_columns = [
    'log R0',
    'log Pi', 'E log Pi', 'e log Pi', 
    'log Lam', 'e log Lam', 'd log Lam',
    'log H', 'log S',
]

dfp.query("Rating >= 4")[select_columns].describe()

dfp.query("Rating == 3")[select_columns].describe()

# So that seems a lot more consistent between the 3 and 4+5-star sub-samples. And it does indeed move everything to larger $\Pi, \Lambda$. 

dfp_by_star = (
    dfp
    .groupby(['Rating'])
    .agg([np.mean, np.std])
)
dfp_by_star[select_columns]

# So the only thing that is still correlated with the star rating is the relative shell thickness (and to lesser extent the breadth).  This *might* be due to errors in `R0`, but I don't think so.

# Write out the peak version to a file for use in later steps. 

dfp.to_csv("mipsgal-arcfit-peak-variants.tab", sep="\t")

# Looking at the distributions, it seems a mild set of cut-offs to eliminate outliers could be:

qmild = "`e log H` < 0.25 and `E log Pi` < 0.4"
dfp.query(f"Rating >= 3 and not ({qmild})").describe()

# So it only eliminates 24, and most of those are 3-star:

dfp.query(f"Rating >= 4 and not ({qmild})")

cols = [
    "log Pi", "log Lam", 
    "log R0", "log H", "log S",
    "E log Pi", 
    "e log Lam", 
    "e log H",
]
g = sns.pairplot(
    dfp.query("Rating >= 3").query(qmild), 
    hue="Rating", 
    vars=cols, 
    palette="magma_r",
    plot_kws=dict(s=50), 
    diag_kind="kde",
)


# We already have the angular thickness column.

dfp[["log R0", "log h", "log H"]].corr()

g = sns.pairplot(
    dfp.query("Rating >= 3").query(qmild), 
    hue="Rating", 
    vars=["log R0", "log h", "log H", "log hin", "log Hin"], 
    palette="magma_r",
    plot_kws=dict(s=20, alpha=0.8, edgecolor="none"), 
    diag_kind="hist",
)
rgrid = np.linspace(0.5, 2.5)
g.axes[1, 0].plot(rgrid, 0.5*rgrid + 0.6)
g.axes[2, 0].plot(rgrid, -0.5*rgrid + 0.7)
g.fig.set_size_inches(15, 12)

# So this is a very strange result: $h$ is positively correlated with $R_0$, with $r^2 \approx 0.7$, and slope of $m \approx 0.5$ (shown by solid line in graph).  Meanwhile, $H$ is negatively correlated with $R_0$, with $r^2 \approx 0.5$ but $m \approx -0.5$.

# Another strange thing is that the correlations are a bit displaced from the Spitzer results.  I need to do another notebook that compares the two. 
# Also, $h$ and $H$ are now weakly negatively correlated. 

dfp["log R0"]


