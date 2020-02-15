# -*- coding: utf-8 -*-
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

# # Comparison of bow shell parameters from WISE and Spitzer
#
# This will include two parts:
#
# 1. Correlation between the same quantity measured for the same source with the two different instruments
# 2. Comparison of the WISE-only sample (generally larger arcs) with the Spitzer sample

import numpy as np
import pandas as pd
from scipy import stats
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_color_codes()
sns.set_context("poster")

pd.set_option('display.precision', 3)
pd.set_option("display.max_columns", 999)
pd.set_option("display.max_rows", 9999)

# ## Set up joint dataframe

# First, read each dataset from its respective file, setting the index to the source ID number (`Seq`).

# +
wdf = pd.read_table(
    "wise-arcfit-peak-variants.tab"
).set_index("Seq")

mdf = pd.read_table(
    "mipsgal-arcfit-peak-variants.tab"
).set_index("Seq")
# -

# Then do an outer join on the index to get a combined dataframe.  The WISE data has suffix ` WISE`, while the Spitzer data has no suffix.

df = wdf.join(
    mdf, 
    how="outer", 
    sort=True, 
    lsuffix=" WISE", 
    rsuffix=""
)

# Replace NaNs with zeros in the rating columns, so we can more easily filter on comparitive ratings.

rating_cols = ["Rating", "Rating WISE"]
df[rating_cols] = df[rating_cols].fillna(value=0)

# Remove columns that we will never need – errors of errors, or columns that are all zero.

to_drop = [c for c in df.columns if "E e" in c or df[c].sum() == 0.0]
df = df.drop(columns=to_drop)


# Define new columns for instrument-corrected widths and breadths.  We subtract the width linearly because the source and PSF both have heavy tails.

# +
def logify(x, dx):
    """Base-10 log of quantity `x` and its error `dx`"""
    return np.log10(x), np.log10(np.e)*dx/x

for suffix, W_ins in ["", 5.5], [" WISE", 12.0]:
    h = 10**df['log h' + suffix]
    # hstar = np.sqrt(h**2 - W_ins**2)
    hstar = h - W_ins
    # An additional error term of 10% of the instrumental width, which we call "E" 
    E_hstar = 0.2*W_ins
    df["log h*" + suffix], df["E log h*" + suffix] = logify(hstar, E_hstar)
    # And retain the original error as "e", but scaled to new value
    df["e log h*" + suffix] = df["e log h" + suffix]*10**(df['log h' + suffix] - df['log h*' + suffix])
    # And do the normalized versions
    df["log H*" + suffix] = df["log h*" + suffix] - df["log R0" + suffix]
    df["e log H*" + suffix] = df["e log h*" + suffix]
    df["E log H*" + suffix] = df["E log h*" + suffix]

    
    s = 10**df['log s' + suffix]
    ss = 10**df['log ss' + suffix]
    # Average of circle and parabola methods
    sstar = 0.5*(s + ss) - W_ins
    # Difference between methods
    d_sstar = 0.5*np.abs(s - ss)
    df["log s*" + suffix], df["d_log s*" + suffix] = logify(sstar, d_sstar)
    # Combine errors from both methods
    df["e log s*" + suffix] = (0.5/1.414)*(df["e log s" + suffix] + df["e log ss" + suffix])
    df["E log s*" + suffix] = (1.0/1.414)*df["E log s" + suffix] # There is no "E log ss"
    # And do the normalized versions
    df["log S*" + suffix] = df["log s*" + suffix] - df["log R0" + suffix]
    df["e log S*" + suffix] = df["e log s*" + suffix]
    df["E log S*" + suffix] = df["E log s*" + suffix]
# -

for suff in "", " WISE":
    df["log K core" + suff], df["e log K core" + suff] = logify(
        df["core_bratio" + suff]/0.5, df["E core_bratio" + suff])
    df["log K wing" + suff], df["e log K wing" + suff] = logify(
        df["wing_bratio" + suff]/1.5, df["E wing_bratio" + suff])
    df["log K wing" + suff] *= -1


# Look at sources where WISE has the better rating

df_wise_best = df.query(
    "`Rating WISE` > Rating and `Rating WISE` >= 3"
)
df_wise_best.describe()

# So there are 114 of them, and all but 15 have no MIPSGAL measurements at all.  Here is a sample of them:

df_wise_best[rating_cols].sample(10).sort_index()

# ## Correlation of $R_0$

# Add a column with the ratio of the radii:

df["log R0 RATIO"] = df["log R0 WISE"] - df["log R0"]

df[["log R0", "log R0 WISE"]].corr()

g = sns.jointplot(
    "log R0", "log R0 RATIO", data=df, 
    kind="scatter", alpha=0.5, edgecolor="none", s=30,
)

# But we really want only those sources that have a good rating in both lists:

df_both_345 = df.query(
    "Rating >= 3 and `Rating WISE` >= 3"
)
df_both_345.describe()

# There are 200 sources in the comparison set. The radius is slightly smaller for WISE: log ratio is $-0.03 \pm 0.07$, which is roughly linear ratio of $0.93 \pm 0.15$.

df_both_345[["log R0", "log R0 WISE"]].corr()

# So, correlation is pretty good.  Here is a plot, with symbol size increasing with WISE star rating:

g = sns.jointplot(
    "log R0", "log R0 RATIO", data=df_both_345, 
    kind="scatter", alpha=0.5, edgecolor="none", 
    s=df_both_345["Rating WISE"]**3,
)
_ = g.ax_joint.axvspan(1.0, 1.3, color="xkcd:red", alpha=0.1, ec="none", zorder=-100)

# The shaded box shows the "mid-sized" sample, with "small-sized" to left and "large-sized" to right.

# It is clear that it is the 3-star sources that have the larger dispersion. 

# Now zoom in on the larger radii. Restrict it to $R_0 > 20''$:

df_both_345_large = df_both_345.query("`log R0` > 1.3")
df_both_345_large.describe()

# This time there are only 40 sources, but the ratio is tighter: $0.97 \pm 0.08$.

df_both_345_large[["log R0", "log R0 WISE"]].corr()

# Correlation is better too, and here is the plot:

g = sns.jointplot(
    "log R0", "log R0 RATIO", data=df_both_345_large, 
    kind="scatter", alpha=0.7, edgecolor="none", 
    s=3*df_both_345_large["Rating WISE"]**3,
)

# So there are 3 outliers with log ratio < -0.09, but apart from that the distribution is very tight.

df_both_345_large.sort_values(
    "log R0 RATIO"
)[["log R0", "log R0 RATIO"]].head(3)

# So, we can remove those:

df_both_345_good = df_both_345_large.query(
    "-0.09 <= `log R0 RATIO` <= 0.09"
)
df_both_345_good.describe()

df_both_345_good[["log R0", "log R0 WISE"]].corr()

# ## Correlations in shape
#
# First an overview.

shape_vars = ["log Pi", "log Lam", "log H*", "log S*"]
shape_vars_wise = [_ + " WISE" for _ in shape_vars]

q_trim_h = "`log h` > 0.5 & `log h WISE` > 0.5 & `log S` > -0.5 & `log H WISE` < 0.0"
g = sns.pairplot(
    df_both_345_good.query(q_trim_h), 
    vars=shape_vars+shape_vars_wise,
    hue="Rating WISE",
    palette="magma_r",
)


# So it all looks good except for the thickness, $h$, which has a few outliers and also seems to be shifted.

# ### Define some library functions for plotting

# +
def total_error(df, var, etypes=["e ", "E ", "d_"]):
    """Combine the e error and the E error in quadrature"""
    err_sqr = np.zeros_like(df[var])
    for e in etypes:
        evar = e + var
        if evar in df:
            err_sqr += df[evar]**2
    return np.sqrt(err_sqr)

    
def get_xy_and_xyerrors(df, var):
    """Return X, dX, Y, dY for the Spitzer (X) and WISE (Y) datasets"""
    varx, vary = var, f"{var} WISE"
    return (df[varx], total_error(df, varx), 
            df[vary], total_error(df, vary))


# -

def plot_with_errors(ax, df, var, v1=None, v2=None, s0=5, escale=0.1, **kwds):
    X, Xe, Y, Ye = get_xy_and_xyerrors(df, var)
    ax.errorbar(X, Y, xerr=Xe, yerr=Ye, ls=" ", alpha=0.1, zorder=-100)
    ax.scatter(X, Y, s=s0/(escale**2 + Xe**2 + Ye**2), 
               c=df["Rating WISE"], vmin=2, vmax=5,
               **{"cmap": "viridis_r", **kwds})
    ax.plot([v1, v2], [v1, v2], c="k", ls="--", alpha=0.5)
    rr = np.corrcoef(X, Y)[0, 1]
    rr_text = f"$r = {rr:.3f}$"
    if "log" in var:
        ratio = 10**(Y - X)
    else:
        ratio = Y/X
    weights = 1.0/np.hypot(Xe, Ye)
    m = np.isfinite(weights)
    ratio_mean = np.average(ratio[m], weights=weights[m])
    ratio_std = np.sqrt(np.average((ratio[m] - ratio_mean)**2, weights=weights[m]))
    ratio_text = fr"WISE/Spitzer = ${ratio_mean:.2f} \pm {ratio_std:.2f}$"
    ax.text(0.98, 0.02, rr_text + "\n" + ratio_text, 
            horizontalalignment='right',
            verticalalignment='bottom', 
            fontsize="small",
            transform=ax.transAxes)
    ax.set(
        xlim=[v1, v2],
        ylim=[v1, v2],
        xlabel="Spitzer " + var,
        ylabel="WISE " + var,
    )   
    ax.set_aspect('equal')


# ### The uncorrected shell thickness and breadth

# This is very strange.  Even for the largest sources, WISE measures a much larger $h$, by a factor of 1.67. And it doesn't get better for the 5-star sources even.

d = df_both_345_good.sort_values(
    "Rating WISE"
).dropna(
    subset=["log h", "log hin", "log h WISE", "log hin WISE"],
)
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_with_errors(ax1, d, "log h", 0.6, 1.8, s0=1, escale=0.05, cmap="twilight_shifted")
plot_with_errors(ax2, d, "log hin", 0.6, 1.8, s0=1, escale=0.05, cmap="twilight_shifted")
# Resolution limits
ax1.axvline(np.log10(5.5), color="k", ls=":")
ax1.axhline(np.log10(12), color="k", ls=":")
ax2.axvline(np.log10(5.5), color="k", ls=":")
ax2.axhline(np.log10(12), color="k", ls=":")
fig.suptitle(
    f"Uncorrected thickness for large-sized overlap sample ($N = {len(d)}$): $R_0 > 20''$",
    fontsize="small",
)
sns.despine(fig)
fig.tight_layout()

# And we even see a deviation in the case of the breadths, which are much larger.

d = df_both_345_good.sort_values(
    "Rating WISE"
).dropna(
    subset=["log s", "log ss", "log s WISE", "log ss WISE"],
)
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_with_errors(ax1, d, "log s", 1.0, 2.0, s0=1, escale=0.05, cmap="twilight_shifted")
plot_with_errors(ax2, d, "log ss", 1.0, 2.0, s0=1, escale=0.05, cmap="twilight_shifted")
# Resolution limits
#ax1.axvline(np.log10(5.5), color="k", ls=":")
#ax1.axhline(np.log10(12), color="k", ls=":")
fig.suptitle(
    f"Uncorrected breadth for large-sized overlap sample ($N = {len(d)}$): $R_0 > 20''$",
    fontsize="small",
)
sns.despine(fig)
fig.tight_layout()

# ### The corrected thickness and breadth

d = df_both_345_good.sort_values(
    "Rating WISE"
).dropna(
    subset=["log h*", "log s*", "log h* WISE", "log s* WISE"],
)
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_with_errors(ax1, d, "log h*", 0.1, 1.4, s0=4, escale=0.05, cmap="twilight_shifted")
plot_with_errors(ax2, d, "log s*", 0.6, 2.2, s0=1, escale=0.05, cmap="twilight_shifted")
# Resolution limits
ax1.axvline(np.log10(5.5), color="k", ls=":")
ax1.axhline(np.log10(12), color="k", ls=":")
ax2.axvline(np.log10(5.5), color="k", ls=":")
ax2.axhline(np.log10(12), color="k", ls=":")
fig.suptitle(
    f"Corrected isophote shape for large-sized overlap sample ($N = {len(d)}$): $R_0 > 20''$",
    fontsize="small",
)
sns.despine(fig)
fig.tight_layout()

d = df_both_345.sort_values("Rating WISE").query(
    "1.0 <= `log R0` <= 1.3"
).dropna(
    subset=["log h*", "log s*", "log h* WISE", "log s* WISE"],
)
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_with_errors(ax1, d, "log h*", 0.1, 1.4, s0=4, escale=0.05, cmap="twilight_shifted")
plot_with_errors(ax2, d, "log s*", 0.6, 2.2, s0=1, escale=0.05, cmap="twilight_shifted")
# Resolution limits
ax1.axvline(np.log10(5.5), color="k", ls=":")
ax1.axhline(np.log10(12), color="k", ls=":")
ax2.axvline(np.log10(5.5), color="k", ls=":")
ax2.axhline(np.log10(12), color="k", ls=":")
sns.despine(fig)
fig.suptitle(
    f"Corrected isophote shape for mid-sized overlap sample ($N = {len(d)}$): $10'' < R_0 < 20''$",
    fontsize="small",
)
fig.tight_layout()

d = df_both_345.sort_values("Rating WISE").query(
    "`log R0` <= 1.0"
).dropna(
    subset=["log h*", "log s*", "log h* WISE", "log s* WISE"],
)
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_with_errors(ax1, d, "log h*", 0.1, 1.4, s0=4.0, escale=0.05, cmap="twilight_shifted")
plot_with_errors(ax2, d, "log s*", 0.6, 2.2, s0=1, escale=0.05, cmap="twilight_shifted")
# Resolution limits
ax1.axvline(np.log10(5.5), color="k", ls=":")
ax1.axhline(np.log10(12), color="k", ls=":")
ax2.axvline(np.log10(5.5), color="k", ls=":")
ax2.axhline(np.log10(12), color="k", ls=":")
sns.despine(fig)
fig.suptitle(
    f"Corrected isophote shape for small-sized overlap sample ($N = {len(d)}$): $R_0 < 10''$",
    fontsize="small",
)
fig.tight_layout()

# ### The planitude and alatude
#
# Conclusions are that:
#
# 1. for the large-R0 sample, both planitude ($r \approx 0.8$) and alatude ($r \approx 0.9$) are well-correlated between the two telescopes
# 2. for mid-sized and small-sized samples, alatude ($r \approx 0.65$) is reasonably well correlated, but planitude ($r \approx 0.3$) is terrible

d = df_both_345_good.sort_values("Rating WISE")
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_with_errors(ax1, d, "log Pi", -0.15, 1.05, s0=1, escale=0.05)
plot_with_errors(ax2, d, "log Lam", -0.15, 1.05, s0=1, escale=0.05)
sns.despine(fig)
fig.suptitle(
    f"Arc shape for large-sized overlap sample ($N = {len(d)}$): $R_0 > 20''$",
    fontsize="small",
)
fig.tight_layout()

d = df_both_345.sort_values("Rating WISE").query("1.0 <= `log R0` <= 1.3")
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_with_errors(ax1, d, "log Pi", -0.15, 1.05, s0=0.5, escale=0.05)
plot_with_errors(ax2, d, "log Lam", -0.15, 1.05, s0=0.5, escale=0.05)
sns.despine(fig)
fig.suptitle(
    f"Arc shape for mid-sized overlap sample ($N = {len(d)}$): $10'' < R_0 < 20''$",
    fontsize="small",
)
fig.tight_layout()

d = df_both_345.sort_values("Rating WISE").query("`log R0` <= 1.0")
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_with_errors(ax1, d, "log Pi", -0.15, 1.05, s0=0.5, escale=0.05)
plot_with_errors(ax2, d, "log Lam", -0.15, 1.05, s0=0.5, escale=0.05)
sns.despine(fig)
fig.suptitle(
    f"Arc shape for small-sized overlap sample ($N = {len(d)}$): $R_0 < 10''$",
    fontsize="small",
)
fig.tight_layout()

# ### The kurtosis
#
# This actually comes out a lot better than I had imagined, especially for the core ratio. 
#
# It would be better to convert them to log values though, and to homogenize the core and wing versions. 
#
# A Gaussian has core = 0.5 and wing = 1.5, with leptokurtic giving small values of core and large values of wing.
#

d = df_both_345_good.sort_values("Rating WISE")
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_with_errors(ax1, d, "core_bratio", 0.0, 1.0, s0=0.05, escale=0.02)
plot_with_errors(ax2, d, "wing_bratio", 1.0, 3.0, s0=0.5, escale=0.05)
sns.despine(fig)
fig.suptitle(
    f"Kurtosis for large-sized overlap sample ($N = {len(d)}$): $R_0 > 20''$",
    fontsize="small",
)
fig.tight_layout()

d = df_both_345_good.sort_values("Rating WISE")
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_with_errors(ax1, d, "log K core", -0.45, 0.35, s0=0.05, escale=0.01)
plot_with_errors(ax2, d, "log K wing", -0.45, 0.35, s0=0.05, escale=0.01)
sns.despine(fig)
fig.suptitle(
    f"Kurtosis for large-sized overlap sample ($N = {len(d)}$): $R_0 > 20''$",
    fontsize="small",
)
fig.tight_layout()

d = df_both_345.sort_values("Rating WISE").query(
    "1.0 <= `log R0` <= 1.3"
).dropna(
    subset=["core_bratio", "wing_bratio", "core_bratio WISE", "wing_bratio WISE"],
)
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_with_errors(ax1, d, "core_bratio", 0.0, 1.0, s0=0.05, escale=0.02)
plot_with_errors(ax2, d, "wing_bratio", 1.0, 3.0, s0=0.5, escale=0.05)
sns.despine(fig)
fig.suptitle(
    f"Kurtosis for mid-sized overlap sample ($N = {len(d)}$): $10'' < R_0 < 20''$",
    fontsize="small",
)
fig.tight_layout()

d = df_both_345.sort_values("Rating WISE").query(
    "1.0 <= `log R0` <= 1.3"
)
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_with_errors(ax1, d, "log K core", -0.45, 0.35, s0=0.05, escale=0.01)
plot_with_errors(ax2, d, "log K wing", -0.45, 0.35, s0=0.05, escale=0.01)
sns.despine(fig)
fig.suptitle(
    f"Kurtosis for medium-sized overlap sample ($N = {len(d)}$): $10'' < R_0 < 20''$",
    fontsize="small",
)
fig.tight_layout()

d = df_both_345.sort_values("Rating WISE").query(
    "`log R0` <= 1.0"
).dropna(
    subset=["core_bratio", "wing_bratio", "core_bratio WISE", "wing_bratio WISE"],
)
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_with_errors(ax1, d, "core_bratio", 0.0, 1.0, s0=0.05, escale=0.02)
plot_with_errors(ax2, d, "wing_bratio", 1.0, 3.0, s0=0.5, escale=0.05)
sns.despine(fig)
fig.suptitle(
    f"Kurtosis for small-sized overlap sample ($N = {len(d)}$): $R_0 < 10''$",
    fontsize="small",
)
fig.tight_layout()

d = df_both_345.sort_values("Rating WISE").query(
    "`log R0` <= 1.0"
)
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_with_errors(ax1, d, "log K core", -0.45, 0.35, s0=0.05, escale=0.01)
plot_with_errors(ax2, d, "log K wing", -0.45, 0.35, s0=0.05, escale=0.01)
sns.despine(fig)
fig.suptitle(
    f"Kurtosis for small-sized overlap sample ($N = {len(d)}$): $10'' < R_0 < 20''$",
    fontsize="small",
)
fig.tight_layout()

# ## Samples of large bow shocks

q_all_large_good = "(Rating >= 3 and `log R0` >= 1.3) or (`Rating WISE` >= 3 and `log R0 WISE` >= 1.3)"
q_all_medium_good = "(Rating >= 3 and 1.0 <= `log R0` <= 1.3) or (`Rating WISE` >= 3 and 1.0 <= `log R0 WISE` <= 1.3)"
q_all_small_good = "(Rating >= 3 and `log R0` <= 1.0) or (`Rating WISE` >= 3 and `log R0 WISE` <= 1.0)"
q_all_large_45 = "(Rating >= 4 and `log R0` >= 1.3) or (`Rating WISE` >= 4 and `log R0 WISE` >= 1.3)"
q_all_medium_45 = "(Rating >= 4 and 1.0 <= `log R0` <= 1.3) or (`Rating WISE` >= 4 and 1.0 <= `log R0 WISE` <= 1.3)"
q_all_small_45 = "(Rating >= 4 and `log R0` <= 1.0) or (`Rating WISE` >= 4 and `log R0 WISE` <= 1.0)"
q_wise_best = "`Rating WISE` > Rating"
mcolumns = [_ for _  in df.columns if not "WISE" in _]
wcolumns = [_ for _  in df.columns if "WISE" in _]

# These samples are all the 3,4,5-star bows with $R_0 > 20''$.  First, the whole lot of them:

df.query(q_all_large_good).describe()

# So, there are 106 in total.  Now, divide into a WISE-best and Spitzer-best sub-samples.  
#
# ### The WISE-best subsample 

df.query(f"({q_all_large_good}) and {q_wise_best}").drop(columns=mcolumns).describe()

# That has $N = 62$.
#
#
# ### The Spitzer-best subsample

df.query(f"({q_all_large_good}) and not {q_wise_best}").drop(columns=wcolumns).describe()


# That has $N = 44$.  

# ### Some more utility functions for graphs

# +
def get_pair_xy_and_xyerrors(df, varx, vary):
    """Return X, dX, Y, dY for two different datasets"""
    return (df[varx], total_error(df, varx), 
            df[vary], total_error(df, vary))




def plot_pair_with_errors(ax, df, suff, varx, vary, v1=None, v2=None, s0=5, escale=0.1, **kwds):
    X, Xe, Y, Ye = get_pair_xy_and_xyerrors(df, varx+suff, vary+suff)
    if len(X):
        ax.errorbar(X, Y, xerr=Xe, yerr=Ye, ls=" ", alpha=0.1, zorder=-100)
        ax.scatter(X, Y, s=s0/(escale**2 + Xe**2 + Ye**2), 
                   c=df["Rating" + suff], vmin=2, vmax=5,
                   **{"cmap": "viridis_r", **kwds})
    ax.plot([v1, v2], [v1, v2], c="k", ls="--", alpha=0.5)
    if "log" in varx:
        ratio = 10**(Y - X)
    else:
        ratio = Y/X
    weights = 1.0/np.hypot(Xe, Ye)
    m = np.isfinite(weights) & np.isfinite(ratio)
    if m.sum():
        rr = np.corrcoef(X[m], Y[m])[0, 1]
        rr_text = f"$r = {rr:.3f}$"
        ratio_mean = np.average(ratio[m], weights=weights[m])
        ratio_std = np.sqrt(np.average((ratio[m] - ratio_mean)**2, weights=weights[m]))
        ratio_text = fr"{vary[4:]}/{varx[4:]} = ${ratio_mean:.2f} \pm {ratio_std:.2f}$"
        ax.text(0.98, 0.02, rr_text + "\n" + ratio_text, 
                horizontalalignment='right',
                verticalalignment='bottom', 
                fontsize="small",
                transform=ax.transAxes)
    ax.set(
        xlim=[v1, v2],
        ylim=[v1, v2],
        xlabel=varx,
        ylabel=vary,
    )   
    ax.set_aspect('equal')


# -

# ### Lambda–Pi for different samples

d_s = df.query(f"({q_all_large_good}) and not {q_wise_best}").sort_values("Rating").drop(columns=wcolumns)
d_w = df.query(f"({q_all_large_good}) and {q_wise_best}").sort_values("Rating WISE").drop(columns=mcolumns)
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_pair_with_errors(ax1, d_s, "", "log Pi", "log Lam", -0.05, 0.85, s0=0.5, escale=0.05)
plot_pair_with_errors(ax2, d_w, " WISE", "log Pi", "log Lam", -0.05, 0.85, s0=0.5, escale=0.05)
ax1.set_title("Spitzer")
ax2.set_title("WISE")
sns.despine(fig)
fig.suptitle(
    f"Arc shape for large-sized samples ($N = {len(d_s)}, {len(d_w)}$): $R_0 > 20''$",
    fontsize="small",
)
fig.tight_layout()

# Look at extreme high values of planitude.

d_s.sort_values("log Pi").tail()

d_w.sort_values("log Pi WISE").tail()

# And low values of planitude

d_s.sort_values("log Pi").head()

# Lowest Spitzer planitude is 587, which is very pointy (v-like) indeed (alatude is not so low).  Next is 014, which has a low-ish S*, possibly indicating R0 too large, but it does look like the right star. Third is 566, which has low alatude too and looks c-like. Next is 221, which is very asymmetrical, should be discarded. Last is 509, which has a large alatude (but uncertain), a smallish S*, and is platykurtic.

d_w.query("`log R0 WISE` >= 1.6").sort_values("log Pi WISE").head()

# Note that I am only selecting $R_0 > 40$ for these, to avoid the poorly resolved ones. The most extreme low-planitude WISE source is 542, which also has a very low dimensionless breadth, $s/R_0$, which makes me think that maybe the star is mis-identified so that $R_0$ is too large. Next one, 363, is indeed pointy. Third one, 054, is more circular.  Fourth one, 371, is point again, possibly with axial bulge. Fifth one, 066, is pointy too.

d_s.sort_values("log Lam").tail()

d_s = df.query(f"({q_all_medium_good}) and not {q_wise_best}").sort_values("Rating")
d_w = df.query(f"({q_all_medium_good}) and {q_wise_best}").sort_values("Rating WISE")
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_pair_with_errors(ax1, d_s, "", "log Pi", "log Lam", -0.05, 0.85, s0=0.5, escale=0.05)
plot_pair_with_errors(ax2, d_w, " WISE", "log Pi", "log Lam", -0.05, 0.85, s0=0.5, escale=0.05)
ax1.set_title("Spitzer")
ax2.set_title("WISE")
sns.despine(fig)
fig.suptitle(
    f"Arc shape for medium-sized samples ($N = {len(d_s)}, {len(d_w)}$): $10'' < R_0 < 20''$",
    fontsize="small",
)
fig.tight_layout()

d_s = df.query(f"({q_all_small_good}) and not {q_wise_best}").sort_values("Rating")
d_w = df.query(f"({q_all_small_good}) and {q_wise_best}").sort_values("Rating WISE")
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_pair_with_errors(ax1, d_s, "", "log Pi", "log Lam", -0.05, 0.85, s0=0.5, escale=0.05)
plot_pair_with_errors(ax2, d_w, " WISE", "log Pi", "log Lam", -0.05, 0.85, s0=0.5, escale=0.05)
ax1.set_title("Spitzer")
ax2.set_title("WISE")
sns.despine(fig)
fig.suptitle(
    f"Arc shape for small-sized samples ($N = {len(d_s)}, {len(d_w)}$): $R_0 < 10''$",
    fontsize="small",
)
fig.tight_layout()

d_s = df.query(f"({q_all_medium_45}) and not {q_wise_best}").sort_values("Rating")
d_w = df.query(f"({q_all_medium_45}) and {q_wise_best}").sort_values("Rating WISE")
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_pair_with_errors(ax1, d_s, "", "log Pi", "log Lam", -0.05, 0.85, s0=0.5, escale=0.05)
plot_pair_with_errors(ax2, d_w, " WISE", "log Pi", "log Lam", -0.05, 0.85, s0=0.5, escale=0.05)
ax1.set_title("Spitzer")
ax2.set_title("WISE")
sns.despine(fig)
fig.suptitle(
    f"Arc shape for medium-sized 4,5-star samples ($N = {len(d_s)}, {len(d_w)}$): $10'' < R_0 < 20''$",
    fontsize="small",
)
fig.tight_layout()

d_s = df.query(f"({q_all_small_45}) and not {q_wise_best}").sort_values("Rating")
d_w = df.query(f"({q_all_small_45}) and {q_wise_best}").sort_values("Rating WISE")
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_pair_with_errors(ax1, d_s, "", "log Pi", "log Lam", -0.05, 0.85, s0=0.5, escale=0.05)
plot_pair_with_errors(ax2, d_w, " WISE", "log Pi", "log Lam", -0.05, 0.85, s0=0.5, escale=0.05)
ax1.set_title("Spitzer")
ax2.set_title("WISE")
sns.despine(fig)
fig.suptitle(
    f"Arc shape for small-sized 4,5-star samples ($N = {len(d_s)}, {len(d_w)}$): $R_0 < 10''$",
    fontsize="small",
)
fig.tight_layout()

# ### Thickness and breadth for different samples

d_s = df.query(f"({q_all_large_good}) and not {q_wise_best}").sort_values("Rating").drop(columns=wcolumns)
d_w = df.query(f"({q_all_large_good}) and {q_wise_best}").sort_values("Rating WISE").drop(columns=mcolumns)
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_pair_with_errors(ax1, d_s, "", "log S*", "log H*", -1.5, 0.8, s0=0.5, escale=0.05)
plot_pair_with_errors(ax2, d_w, " WISE", "log S*", "log H*", -1.5, 0.8, s0=0.5, escale=0.05)
ax1.set_title("Spitzer")
ax2.set_title("WISE")
sns.despine(fig)
fig.suptitle(
    f"Breadth and thickness for large-sized samples ($N = {len(d_s)}, {len(d_w)}$): $R_0 > 20''$",
    fontsize="small",
)
fig.tight_layout()

d_s.dropna(subset=["log H*"]).sort_values("log H*").tail()

d_w.dropna(subset=["log H* WISE"]).sort_values("log H* WISE").tail()



d_s = df.query(f"({q_all_medium_45}) and not {q_wise_best}").sort_values("Rating")
d_w = df.query(f"({q_all_medium_45}) and {q_wise_best}").sort_values("Rating WISE")
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_pair_with_errors(ax1, d_s, "", "log S*", "log H*", -1.5, 0.8, s0=0.5, escale=0.05)
plot_pair_with_errors(ax2, d_w, " WISE", "log S*", "log H*", -1.5, 0.8, s0=0.5, escale=0.05)
ax1.set_title("Spitzer")
ax2.set_title("WISE")
sns.despine(fig)
fig.suptitle(
    f"Breadth and thickness for medium-sized samples ($N = {len(d_s)}, {len(d_w)}$): $10'' < R_0 < 20''$",
    fontsize="small",
)
fig.tight_layout()

d_s = df.query(f"({q_all_small_45}) and not {q_wise_best}").sort_values("Rating")
d_w = df.query(f"({q_all_small_45}) and {q_wise_best}").sort_values("Rating WISE")
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7))
plot_pair_with_errors(ax1, d_s, "", "log S*", "log H*", -1.5, 0.8, s0=0.5, escale=0.05)
plot_pair_with_errors(ax2, d_w, " WISE", "log S*", "log H*", -1.5, 0.8, s0=0.5, escale=0.05)
ax1.set_title("Spitzer")
ax2.set_title("WISE")
sns.despine(fig)
fig.suptitle(
    f"Breadth and thickness for small-sized samples ($N = {len(d_s)}, {len(d_w)}$): $R_0 < 10''$",
    fontsize="small",
)
fig.tight_layout()

# ## Look into the thickness problem

selected_vars = ["log R0", "Rating", "Rating WISE", "log hin", "log hin WISE", "log h", "log h WISE", "log h*", "log h* WISE"]
df_both_345_good[selected_vars].sort_values("log R0")


