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

# And plot the distribution of $\Pi$ and $\Lambda$ with error bars. 

# %matplotlib inline
from matplotlib import pyplot as plt
import seaborn as sns

fig, ax = plt.subplots(figsize=(10, 6))
ax.errorbar(x, y, xerr=xe, yerr=ye, ls=" ", alpha=0.5)
ax.set(
    xlim=[0, None], ylim=[0, None],
    xlabel=r"$\Pi$", ylabel=r"$\Lambda$",
)
ax.set_aspect("equal")

lm = linmix.LinMix(x, y, xe, ye, K=2)

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

pearsonr(x, y)

# So this is for the 4,5-star sources.  We get a posterior distribution for the correlation coefficient of $r = 0.57 \pm 0.1$, assuming this is what `lm.chain["corr"]` is.  This is significantly larger than the correlation in observed variables: $r = 0.44$ (see previous cell).

tab[m][useful_cols].to_pandas().describe()

# So note that the observed variable $\Pi$ (`Rc`) is $2.026 \pm 0.777$, whereas the latent variable $\xi$ from the chain posterior is $2.022 \pm 0.770$, which is hardly any different.  
#
# This is probably because the observaional error bars on $\Pi$ are so small, $0.089 \pm 0.049$.  We should probably increase them, by consideration of different $\theta$ ranges for the circle fits.  In Paper 0 we found values of 0.13 to 0.38, so 0.25 might be a better bet. We could multiply `Rc_sigma` by 3 and see what happens

lm2 = linmix.LinMix(x, y, 3*xe, ye, K=2, seed=443)

lm2.run_mcmc()

dfchain2 = pd.DataFrame.from_records(
    lm2.chain.tolist(), 
    columns=lm2.chain.dtype.names
)
dfchain2.describe()

# So now we get $2.015 \pm 0.718$, which is hardly any different from before, although the width is very slightly smaller $0.718 \pm 0.084$ instead of $0.771 \pm 0.086$.
#
# The uncertainty on the slope has increased, as one would hope, from $0.157$ to $0.167 \pm 0.04$. 
#


