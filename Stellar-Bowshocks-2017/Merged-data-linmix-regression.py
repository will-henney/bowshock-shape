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

import numpy as np
import linmix

# Read the data of the bow shock shapes from the multiple fit variants

from astropy.table import Table
combo_file = 'mipsgal-arcfit-peak-variants.tab'
tab = Table.read(combo_file, format='ascii.tab')
tab

# Combine the internal and external errors

tab['log Lam sigma'] = np.hypot(tab['e log Lam'], tab['E log Lam'])
tab['log Pi sigma'] = np.hypot(tab['e log Pi'], tab['E log Pi'])

# Restrict to only the 4 and 5 star sources.  Note that I originally tried this with all the sources, but linmix ran all night and failed to produce any output.

m = tab["Rating"] >= 4

useful_cols = ['log Pi', 'log Pi sigma', 'log Lam', 'log Lam sigma']
tab[m][useful_cols]

# So that is 144 points.  We extract the columns from the table as individual variables to use with linmix. 

X, Xe, Y, Ye = [tab[m][_] for _ in ['log Pi', 'log Pi sigma', 'log Lam', 'log Lam sigma']]

# And plot the distribution of $\Pi$ and $\Lambda$ with error bars. 

# %matplotlib inline
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_context("talk")

fig, ax = plt.subplots(figsize=(10, 10))
ax.errorbar(X, Y, xerr=Xe, yerr=Ye, ls=" ", alpha=0.3)
ax.set(
    xlim=[-0.2, 0.8], ylim=[-0.2, 0.8],
    xlabel=r"$\log_{10}\, \Pi$", ylabel=r"$\log_{10}\, \Lambda$",
)
ax.set_aspect("equal")
sns.despine()

# And do the MCMC regression

lm = linmix.LinMix(X, Y, Xe, Ye, K=2)

lm.run_mcmc()

# The result from the regression is given as `lm.chain`, which consiste of 5000 sampls ver the posterior distribution of the fit parameters.  This is a recarray, which we convert to a pandas data table for convenience. (Note the `.tolist()` trick, which is necessary since some of the columns are 2-vectors. 

import pandas as pd

dfchain = pd.DataFrame.from_records(
    lm.chain.tolist(), 
    columns=lm.chain.dtype.names
)
dfchain

# The intrinsic relation between the latent variables $\xi, \eta$ is 
# $$
# \eta = \alpha + \beta \xi \pm \sigma
# $$
# with paremeters given in the first three columns:
# 1. `alpha` is the intercept of the regression line
# 2. `beta` is the slope of the regression line
# 3. `sigsqr` is the intrinsic variance $\sigma^2$ of $\eta$ around the regression line.
#
# Then there is the derived Gaussian mixture intrinsic distribution of $\xi$ (we asked for two components by giving `K=2`). 
#
# 4. `pi` is height of component
# 5. `mu` is location (central $\xi$ value) of component
# 6. `tausqr` is square of rms width of component
#
# Each of these parameters is a 2-tuple for the two Gaussian mixture components. 

# Then, the parameters $\mu_0$, $u^2$, $w^2$ (columns 7 `mu0`, 8 `usqr`, 9 `wsqr`) are another level of hyper-parameters that are used in constructing the priors for the Gaussian mixture (hyperprior). 

# Finally, we have the mean and std of $\xi$ (columns 10 `ximean`, 11 `xisig`) and the intrinsic correlation of the latent variables $r(\eta, \xi)$ (column 12 `corr`).

dfchain.describe()

sns.pairplot(dfchain.query("xisig < 0.2"), plot_kws=dict(alpha=0.1, s=2), 
             vars=["alpha", "beta", "corr", "ximean", "xisig"])

from scipy.stats import pearsonr

pearsonr(X, Y)

# So this is for the 4,5-star sources in log space.  We get a posterior distribution for the correlation coefficient of $r = 0.91 \pm 0.07$.  This is significantly larger than the correlation in observed log variables: $r = 0.61$ (see previous cell).

pd.DataFrame({"X": X, "Xe": Xe, "Y": Y, "Ye": Ye}).describe()

# So note that the observed variable $\log_{10} \Pi$ (`X`) is $0.348 \pm 0.188$, whereas the latent variable $\xi$ from the chain posterior is $0.317 \pm 0.139$, which is narrower and very slightly shifted.  
#
# The value of $\sigma$ is $0.027 \pm 0.010$, which is the intrinsic spread in the latent $\eta(\xi)$ (or true $\log\Lambda(\log\Pi)$) relation. This is small compared with $\beta \sigma_\xi = 0.06$, which is why the correlation is so high. 
#
# Unlike in the previous notebook, we don't need to do a run with larger X errors since we already have a good estimate of the errors. 

# Plot the regression (note that the biggest points have the smallest errors). 

# +
vmin, vmax = -0.2, 1.0
xgrid = np.linspace(vmin, vmax, 200)


fig, ax = plt.subplots(figsize=(10, 10))

ax.errorbar(X, Y, xerr=Xe, yerr=Ye, ls=" ", elinewidth=0.4, alpha=1.0, c="k")
ax.scatter(X, Y, marker=".", s=20/np.hypot(Xe, Ye))
# The original fit
ax.plot(xgrid, dfchain["alpha"].mean() + xgrid*dfchain["beta"].mean(), 
        '-', c="k")
for samp in lm.chain[::20]:
    ax.plot(xgrid, samp["alpha"] + xgrid*samp["beta"], 
        '-', c="r", alpha=0.2, lw=0.1)

ax.set(
    xlim=[vmin, vmax], ylim=[vmin, vmax],
    xlabel=r"$\log_{10}\, \Pi$", ylabel=r"$\log_{10}\, \Lambda$",
)
ax.set_aspect("equal")
sns.despine()
# -

# So this is now consistent with 0.5, as expected for parabolas, without artificial inflation of the error.  It also has a less uncertainty in the slope than last time.

# Now try the 3-star sources

m3 = tab["Rating"] == 3
X3, Xe3, Y3, Ye3 = [tab[m3][_] for _ in ['log Pi', 'log Pi sigma', 'log Lam', 'log Lam sigma']]

mgood = np.isfinite(Xe3) & np.isfinite(Ye3)
X3 = X3[mgood]
Y3 = Y3[mgood]
Xe3 = Xe3[mgood]
Ye3 = Ye3[mgood]

fig, ax = plt.subplots(figsize=(10, 10))
ax.errorbar(X3, Y3, xerr=Xe3, yerr=Ye3, ls=" ", c="r", alpha=0.5)
ax.set(
    xlim=[-0.25, 1.05], ylim=[-0.25, 1.05],
    xlabel=r"$\log_{10}\, \Pi$", ylabel=r"$\log_{10}\, \Lambda$",
)
ax.set_aspect("equal")
sns.despine()

lm3 = linmix.LinMix(X3, Y3, Xe3, Ye3, K=2)

lm3.run_mcmc()

dfchain3 = pd.DataFrame.from_records(
    lm3.chain.tolist(), 
    columns=lm3.chain.dtype.names
)
dfchain3.describe()

pd.DataFrame({"X": X3, "Xe": Xe3, "Y": Y3, "Ye": Ye3}).describe()

pearsonr(X3, Y3)

# So, again we get a nice slope, albeit a little lower than the 4,5-star. And the correlation is consistent within 1 sigma. 

# +
vmin, vmax = -0.2, 1.1
xgrid = np.linspace(vmin, vmax, 200)


fig, ax = plt.subplots(figsize=(10, 10))

# The 3-star points
ax.errorbar(X3, Y3, xerr=Xe3, yerr=Ye3, ls=" ", elinewidth=0.4, alpha=0.3, c="k")
ax.scatter(X3, Y3, marker=".", c="c", s=20/np.hypot(Xe, Ye), alpha=0.5)

# The 4,5-star points
ax.errorbar(X, Y, xerr=Xe, yerr=Ye, ls=" ", elinewidth=0.4, alpha=0.3, c="k")
ax.scatter(X, Y, marker=".", c="y", s=20/np.hypot(Xe, Ye))

# The original fit
ax.plot(xgrid, dfchain3["alpha"].mean() + xgrid*dfchain3["beta"].mean(), 
        '-', c="k")
for samp in lm3.chain[::20]:
    ax.plot(xgrid, samp["alpha"] + xgrid*samp["beta"], 
        '-', c="r", alpha=0.2, lw=0.1)
# The 4,5-star fit
ax.plot(xgrid, dfchain["alpha"].mean() + xgrid*dfchain["beta"].mean(), 
        '-', c="k")
for samp in lm.chain[::20]:
    ax.plot(xgrid, samp["alpha"] + xgrid*samp["beta"], 
        '-', c="m", alpha=0.2, lw=0.1)


ax.set(
    xlim=[vmin, vmax], ylim=[vmin, vmax],
    xlabel=r"$\log_{10}\, \Pi$", ylabel=r"$\log_{10}\, \Lambda$",
)
ax.set_aspect("equal")
sns.despine()
# -

dfchain3p45 = pd.concat([dfchain, dfchain3], keys=["4,5-star", "3-star"])

dfchain3p45

# Using `pd.concat` on the two dataframes gives a multi-level index.  In order to get seaborn's pairplot to distinguish the two, we need to switch that to a categorical column, which we can do with `.reset_index()`: 

dfchain3p45.reset_index(level=0)

dfchain3p45["sigma"] = np.sqrt(dfchain3p45["sigsqr"])

g = sns.pairplot(dfchain3p45[0::50].reset_index().query("xisig < 0.3 and 0.0 < sigsqr < 0.004"), 
             #kind="reg", 
             diag_kind="hist",
             markers=".",
             plot_kws=dict(marker=".", ec="none", alpha=0.5, s=100), 
             vars=["alpha", "beta", "sigma", "corr", "ximean", "xisig"], 
             hue="level_0",
             hue_order=["3-star", "4,5-star",]
            )
g.axes[2, 0].set(ylim=[0.0, 0.06])

# So the above shows that the results from the two data partitions are *now completely different*.  All except for the mean of $\xi$ (the latent planitude)
#
#
# The 95% credibility levels on $\beta$ and the other parameters are:
#

dfchain.quantile([0.025, 0.5, 0.975])

dfchain3.quantile([0.025, 0.5, 0.975])

# Or 90% confidence levels of

dfchain.quantile([0.05, 0.95])

# Note that the parabola family should have $\Lambda = (2 \Pi)^{1/2}$, so $\alpha, \beta = 0.15, 0.5$.  The slope is consistent with the 4+5-star data, but the intercept is excluded at the 95% level.  The 3-star data gets the slope, but not the intercept.

# ## Application to isophotal shape variables

# What we have at the moment are thickness $h$ and breadth $d\theta$ (which we want to replace with a radius eventually). 

# ### Planitude versus angular breadth

# Start off with $\Pi$ and $d\theta$.  We won't take log of the breadth, for ease of comparison with what I have already done. 

# +
# Add extra 10 deg uncertainty
tab['breadth50 sigma'] = np.hypot(tab['d_breadth50'], 10.0)

# Use 3 to 5 star
m = tab["Rating"] >= 3

X, Xe, Y, Ye = [tab[m][_] for _ in ['log Pi', 'log Pi sigma', 'breadth50', 'breadth50 sigma']]
# -

mgood = np.isfinite(Xe) & np.isfinite(Ye)
X = X[mgood]
Y = Y[mgood]
Xe = Xe[mgood]
Ye = Ye[mgood]

fig, ax = plt.subplots(figsize=(10, 10))
ax.errorbar(X, Y, xerr=Xe, yerr=Ye, ls=" ", c="m", alpha=0.2)
ax.set(
    xlim=[-0.25, 1.05], ylim=[0, 150],
    xlabel=r"$\log_{10}\, \Pi$", ylabel=r"$d\theta$, degrees",
)
sns.despine()

lm_Rc_thb = linmix.LinMix(X, Y, Xe, Ye, K=2)
lm_Rc_thb.run_mcmc()

dfchain_Rc_thb = pd.DataFrame.from_records(
    lm_Rc_thb.chain.tolist(), 
    columns=lm_Rc_thb.chain.dtype.names
)
dfchain_Rc_thb.describe()

dfchain_Rc_thb.head()

# +
vmin, vmax = -0.25, 1.05
xgrid = np.linspace(vmin, vmax, 200)


fig, ax = plt.subplots(figsize=(10, 10))

ax.errorbar(X, Y, xerr=Xe, yerr=Ye, ls=" ", elinewidth=0.4, alpha=1.0, c="k")
ax.scatter(X, Y, marker=".", c="g", s=20/Xe)
# The original fit
ax.plot(xgrid, dfchain_Rc_thb["alpha"].mean() + xgrid*dfchain_Rc_thb["beta"].mean(), 
        '-', c="k")
for samp in lm_Rc_thb.chain[::20]:
    ax.plot(xgrid, samp["alpha"] + xgrid*samp["beta"], 
        '-', c="r", alpha=0.2, lw=0.1)

ax.set(
    xlim=[vmin, vmax], ylim=[0.0, 150.0],
    xlabel=r"$\log_{10}\, \Pi$", ylabel=r"$d\theta$, degrees",
)
sns.despine()
# -

pearsonr(X, Y)

# For some reason, the correlation in observables is less than before: -0.20 instead of -0.47. The latent correlation is $-0.45 \pm 0.15$, which is not back to being significant again.

# Here are the 95% confidence limits:

dfchain_Rc_thb.quantile([0.025, 0.975])

# And 99% confidence

dfchain_Rc_thb.quantile([0.005, 0.995])

# And 99.9% for good measure:

dfchain_Rc_thb.quantile([0.0005, 0.9995])

# So, a latent correlation of zero is excluded at the 99% level, but not at the 99.9% level. 

from scipy.stats import percentileofscore

percentileofscore(dfchain_Rc_thb["corr"], 0.0)

# So, it is at the 99.87 level, or $p = 0.0013$, which is pretty good actually.


