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

# Read the data of the bow shock shapes

from astropy.table import Table
combo_file = 'mipsgal-arcfit.tab'
tab = Table.read(combo_file, format='ascii.tab')
tab

for col in ['Rc', 'Rc_sigma', 'R90p', 'R90p_sigma', 'R90n', 'R90n_sigma']:
      tab[col] /= tab['R0_fit']

# For the errors in $\Lambda$ and $\Pi$ we need to use the uncertainties in R0 too

tab['R90'] = 0.5*(tab['R90p'] + tab['R90n'])
tab['R90_sigma'] = np.sqrt((tab['R0_sigma']/tab['R0_fit'])**2 + 0.5*(tab['R90n_sigma']**2 + tab['R90p_sigma']**2) )
tab['R90_asym'] =  0.5*np.abs(tab['R90p'] - tab['R90n'])
tab['Rc_sigma'] = np.sqrt((tab['R0_sigma']/tab['R0_fit'])**2 + tab['Rc_sigma']**2)

useful_cols = ['R0_fit', 'Rc', 'Rc_sigma', 'R90', 'R90_sigma', 'R90_asym']
for col in useful_cols:
    tab[col] = np.round(tab[col], 3)

# Restrict to only the 4 and 5 star sources.  Note that I originally tried this with all the sources, but linmix ran all night and failed to produce any output.

m = tab["Rating"] >= 4

tab[m][useful_cols]

# So that is 94 points.  We extract the columns from the table as indivisual variables to use with linmix. 

t = tab[m][['Rc', 'Rc_sigma', 'R90', 'R90_sigma']]
x, xe, y, ye = [np.array(_) for _ in zip(*t.as_array().data)]

# Convert to a log scale.  Calculate the errors explicitly, using centered difference so they are symmetric in log  space.  Use capital letters for log quantities.

X, Y = np.log10(x), np.log10(y)
Xe = 0.5*(np.log10(x + xe) - np.log10(x - xe))
Ye = 0.5*(np.log10(y + ye) - np.log10(y - ye))

# And plot the distribution of $\Pi$ and $\Lambda$ with error bars. 

# %matplotlib inline
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_context("talk")

fig, ax = plt.subplots(figsize=(10, 10))
ax.errorbar(X, Y, xerr=Xe, yerr=Ye, ls=" ", alpha=0.5)
ax.set(
    xlim=[-0.5, 0.8], ylim=[-0.5, 0.8],
    xlabel=r"$\log_{10}\, \Pi$", ylabel=r"$\log_{10}\, \Lambda$",
)
ax.set_aspect("equal")
sns.despine()

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

sns.pairplot(dfchain, plot_kws=dict(alpha=0.1, s=2), 
             vars=["alpha", "beta", "corr"])

from scipy.stats import pearsonr

pearsonr(X, Y)

# So this is for the 4,5-star sources in log space.  We get a posterior distribution for the correlation coefficient of $r = 0.66 \pm 0.09$.  This is significantly larger than the correlation in observed log variables: $r = 0.53$ (see previous cell).

pd.DataFrame({"X": X, "Xe": Xe, "Y": Y, "Ye": Ye}).describe()

# So note that the observed variable $\log_{10} \Pi$ (`X`) is $0.279 \pm 0.151$, whereas the latent variable $\xi$ from the chain posterior is $0.283 \pm 0.149$, which is hardly any different.  
#
# This is probably because the observational error bars on $\Pi$ are so small, $0.021 \pm 0.016$.  We should probably increase them, by consideration of different $\theta$ ranges for the circle fits.  In Paper 0 we found values of 0.13 to 0.38 in linear space, so 0.25 might be a better bet, which corresponds to 0.1 in log space. We could add 0.1 to all the `Xe` and see what happens

lm2 = linmix.LinMix(X, Y, Xe+0.1, Ye, K=2, seed=443)

lm2.run_mcmc()

dfchain2 = pd.DataFrame.from_records(
    lm2.chain.tolist(), 
    columns=lm2.chain.dtype.names
)
dfchain2.describe()

# So now we get $0.290 \pm 0.101$ for the `xi` distribution, which is 30% narrower than before. And we get a much higher latent correlation of $r = 0.89 \pm 0.10$, corresponding to $\sigma = 0.026$ (cf the observed dispersion in `Y` of 0.07). 
#
# Also, the slope $\beta$  has changed significantly from $0.25 \pm 0.04$ to $0.53 \pm 0.12$.  

# Plot the regression

# +
vmin, vmax = -0.5, 0.8
xgrid = np.linspace(vmin, vmax, 200)


fig, ax = plt.subplots(figsize=(10, 10))

ax.errorbar(X, Y, xerr=Xe+0.1, yerr=Ye, ls=" ", elinewidth=0.2, alpha=1.0, c="k")
ax.errorbar(X, Y, xerr=Xe, yerr=Ye, ls=" ", elinewidth=0.4, alpha=1.0, c="k")
ax.scatter(X, Y, marker=".")
# The original fit
ax.plot(xgrid, dfchain["alpha"].mean() + xgrid*dfchain["beta"].mean(), 
        '-', c="k")
for samp in lm.chain[::20]:
    ax.plot(xgrid, samp["alpha"] + xgrid*samp["beta"], 
        '-', c="r", alpha=0.2, lw=0.1)
# The fit with increased X errors
ax.plot(xgrid, dfchain2["alpha"].mean() + xgrid*dfchain2["beta"].mean(), 
        '-', c="k")
for samp in lm2.chain[::20]:
    ax.plot(xgrid, samp["alpha"] + xgrid*samp["beta"], 
        '-', c="m", alpha=0.2, lw=0.1)

ax.set(
    xlim=[vmin, vmax], ylim=[vmin, vmax],
    xlabel=r"$\log_{10}\, \Pi$", ylabel=r"$\log_{10}\, \Lambda$",
)
ax.set_aspect("equal")
sns.despine()
# -

# So that graphically shows the difference between the two fits.  The fit with just the nominal $\Pi$ errors (red) is far too shallow. Allowing for an extra 0.1 dex uncertainty in $\Pi$ (magenta) gives a much steeper slope, which is now consistent with 0.5, as expected for parabolas.  It also has a brader uncertainty in the slope, possible even bimodal.

# The worrying thing is that there is nothing in the method that warns us that the original fit was bad.  I suppose there must be a $\chi^2$ that we could have calculated.  

# Now look at the posterior of fit parameters for the new version.

sns.pairplot(dfchain2, plot_kws=dict(alpha=0.01, s=1), 
             vars=["alpha", "beta", "corr"])

# Now try the 3-star sources

m3 = tab["Rating"] == 3
t = tab[m3][['Rc', 'Rc_sigma', 'R90', 'R90_sigma']]
x, xe, y, ye = [np.array(_) for _ in zip(*t.as_array().data)]
X3, Y3 = np.log10(x), np.log10(y)
Xe3 = 0.5*(np.log10(x + xe) - np.log10(x - xe))
Ye3 = 0.5*(np.log10(y + ye) - np.log10(y - ye))

mgood = np.isfinite(Xe3) & np.isfinite(Ye3)
X3 = X3[mgood]
Y3 = Y3[mgood]
Xe3 = Xe3[mgood]
Ye3 = Ye3[mgood]

fig, ax = plt.subplots(figsize=(10, 10))
ax.errorbar(X3, Y3, xerr=Xe3, yerr=Ye3, ls=" ", c="r", alpha=0.5)
ax.set(
    xlim=[-0.5, 0.8], ylim=[-0.5, 0.8],
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

# So, again we get a flat slope. Let us try with increasing the `Xe3` errors by 0.1 again 

lm3e = linmix.LinMix(X3, Y3, Xe3+0.1, Ye3, K=2)
lm3e.run_mcmc()

dfchain3e = pd.DataFrame.from_records(
    lm3e.chain.tolist(), 
    columns=lm3e.chain.dtype.names
)
dfchain3e.describe()

# That is good.  So with this sample too, we now recover a slope of 0.5

# +
vmin, vmax = -0.5, 0.8
xgrid = np.linspace(vmin, vmax, 200)


fig, ax = plt.subplots(figsize=(10, 10))

ax.errorbar(X3, Y3, xerr=Xe3+0.1, yerr=Ye3, ls=" ", elinewidth=0.2, alpha=1.0, c="k")
ax.errorbar(X3, Y3, xerr=Xe3, yerr=Ye3, ls=" ", elinewidth=0.4, alpha=1.0, c="k")
ax.scatter(X3, Y3, marker=".", c="g")
# The original fit
ax.plot(xgrid, dfchain3["alpha"].mean() + xgrid*dfchain3["beta"].mean(), 
        '-', c="k")
for samp in lm3.chain[::20]:
    ax.plot(xgrid, samp["alpha"] + xgrid*samp["beta"], 
        '-', c="r", alpha=0.2, lw=0.1)
# The fit with increased X errors
ax.plot(xgrid, dfchain3e["alpha"].mean() + xgrid*dfchain3e["beta"].mean(), 
        '-', c="k")
for samp in lm3e.chain[::20]:
    ax.plot(xgrid, samp["alpha"] + xgrid*samp["beta"], 
        '-', c="m", alpha=0.2, lw=0.1)

ax.set(
    xlim=[vmin, vmax], ylim=[vmin, vmax],
    xlabel=r"$\log_{10}\, \Pi$", ylabel=r"$\log_{10}\, \Lambda$",
)
ax.set_aspect("equal")
sns.despine()
# -

dfchain3p45 = pd.concat([dfchain2, dfchain3e], keys=["4,5-star", "3-star"])

dfchain3p45

# Using `pd.concat` on the two dataframes gives a multi-level index.  In order to get seaborn's pairplot to distinguish the two, we need to switch that to a categorical column, which we can do with `.reset_index()`: 

dfchain3p45.reset_index(level=0)

sns.pairplot(dfchain3p45[::30].reset_index(), 
             #kind="reg", 
             diag_kind="hist",
             markers=".",
             plot_kws=dict(marker=".", ec="none", alpha=0.3, s=100), 
             vars=["alpha", "beta", "corr", "ximean", "xisig"], 
             hue="level_0",
             hue_order=["3-star", "4,5-star",]
            )

# So the above shows that the results from the two data partitions are very similar.  All except for the mean of $\xi$ (the latent planitude)
#
#
# The 95% confidence levels on $\beta$ and the other parameters are:
#

dfchain2.quantile([0.025, 0.975])

dfchain3e.quantile([0.025, 0.975])

# Or 90% confidence levels of

dfchain2.quantile([0.05, 0.95])

# ## Application to isophotal shape variables

# What we have at the moment are thicknes $h$ and breadth $d\theta$ (which we want to replace with a radius eventually). 

# Start off with $\Pi$ and $d\theta$.  We won't take log of the breadth, for ease of comparison with what I have already done. 

# +
# Half-width angular breadth
tab["thb"] = 0.5*(tab["th_p"] - tab["th_m"])
# Use the asymmetry as estimate of the uncertainty
tab["thb_sigma"] = 0.5*(tab["th_p"] + tab["th_m"])

# Use 3 to 5 star
m = tab["Rating"] >= 3

# Convert to standard form (only log10 on x axis)
t = tab[m][['Rc', 'Rc_sigma', 'thb', 'thb_sigma']]
x, xe, y, ye = [np.array(_) for _ in zip(*t.as_array().data)]
X, Y = np.log10(x), y
Xe = 0.5*(np.log10(x + xe) - np.log10(x - xe))
Ye = ye

# -

mgood = np.isfinite(Xe) & np.isfinite(Ye)
X = X[mgood]
Y = Y[mgood]
Xe = Xe[mgood]
Ye = Ye[mgood]

fig, ax = plt.subplots(figsize=(10, 10))
ax.errorbar(X, Y, xerr=Xe, yerr=Ye, ls=" ", c="m", alpha=0.5)
ax.set(
    xlim=[-0.5, 0.8], ylim=[0, 150],
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

lm_Rc_thb_e = linmix.LinMix(X, Y, Xe + 0.1, Ye, K=2)
lm_Rc_thb_e.run_mcmc()

dfchain_Rc_thb_e = pd.DataFrame.from_records(
    lm_Rc_thb_e.chain.tolist(), 
    columns=lm_Rc_thb_e.chain.dtype.names
)
dfchain_Rc_thb_e.describe()

dfchain_Rc_thb_e.head()

# +
vmin, vmax = -0.5, 0.8
xgrid = np.linspace(vmin, vmax, 200)


fig, ax = plt.subplots(figsize=(10, 10))

ax.errorbar(X, Y, xerr=Xe+0.1, yerr=Ye, ls=" ", elinewidth=0.2, alpha=1.0, c="k")
ax.errorbar(X, Y, xerr=Xe, yerr=Ye, ls=" ", elinewidth=0.4, alpha=1.0, c="k")
ax.scatter(X, Y, marker=".", c="g")
# The original fit
ax.plot(xgrid, dfchain_Rc_thb["alpha"].mean() + xgrid*dfchain_Rc_thb["beta"].mean(), 
        '-', c="k")
for samp in lm_Rc_thb.chain[::20]:
    ax.plot(xgrid, samp["alpha"] + xgrid*samp["beta"], 
        '-', c="r", alpha=0.2, lw=0.1)
# The fit with increased X errors
ax.plot(xgrid, dfchain_Rc_thb_e["alpha"].mean() + xgrid*dfchain_Rc_thb_e["beta"].mean(), 
        '-', c="k")
for samp in lm_Rc_thb_e.chain[::20]:
    ax.plot(xgrid, samp["alpha"] + xgrid*samp["beta"], 
        '-', c="m", alpha=0.2, lw=0.1)

ax.set(
    xlim=[vmin, vmax], ylim=[0.0, 150.0],
    xlabel=r"$\log_{10}\, \Pi$", ylabel=r"$d\theta$, degrees",
)
sns.despine()
# -

pearsonr(X, Y)

# For some reason, the correlation in observables is less than before: -0.34 instead of -0.47. The latent correlation is $-0.37 \pm 0.15$, which is not very well constrained.  

dfchain_Rc_thb_both = pd.concat([dfchain_Rc_thb, dfchain_Rc_thb_e], keys=["small sigma", "large sigma"])

sns.pairplot(dfchain_Rc_thb_both[::30].reset_index(), 
             #kind="reg", 
             diag_kind="hist",
             markers=".",
             plot_kws=dict(marker=".", ec="none", alpha=0.3, s=100), 
             vars=["alpha", "beta", "corr", "ximean", "xisig"], 
             hue="level_0",
             hue_order=["large sigma", "small sigma",]
            )

# The above shows graphically the large difference between the original (large sigma, orange points) and amplified (small sigma, blue points) observational errors on $\Pi$.  With large sigma, the posterior distributions of $\alpha, \beta, r$ are broadened, while $\sigma_\xi$ is shifted to a smaller value. Interestingly, even $\langle \xi \rangle$ moves a bit. 
#
# Here are the 95% confidence limits:
#
#

dfchain_Rc_thb.quantile([0.025, 0.975])

dfchain_Rc_thb_e.quantile([0.025, 0.975])

# And 99% confidence for the large sigma version:

dfchain_Rc_thb_e.quantile([0.005, 0.995])

# And 99.9% for good measure:

dfchain_Rc_thb_e.quantile([0.0005, 0.9995])

# So, a latent correlation of zero is excluded at the 99% level, but not at the 99.9% level. 

from scipy.stats import percentileofscore

percentileofscore(dfchain_Rc_thb_e["corr"], 0.0)

# So, it is at the 99.81 level, or $p = 0.0019$

percentileofscore(dfchain_Rc_thb["corr"], 0.0)

# Whereas for the original sigmas, none of the posterior samples had a positive correlation.  So $p < 0.0001$ since we had 10000 samples.

len(dfchain_Rc_thb)


