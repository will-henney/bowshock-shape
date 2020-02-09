import numpy as np
from scipy.stats import linregress, pearsonr
from sklearn.linear_model import TheilSenRegressor 
from matplotlib import pyplot as plt


def robust_linregress(x, y):
    """Simple interface to other regressors

    This has same signature as `scipy.stats.linregress`
    """

    # The Pearson correlation is independent of the regression
    r, p = pearsonr(x, y)

    # Fit a linear regression y = a x + b
    reg = TheilSenRegressor(random_state=443)
    X = x[:, None]
    reg.fit(X, y)
    a, = reg.coef_
    b = reg.intercept_
    da = 0.0
    return a, b, r, p, da


def plot_regression(ax, xdata, ydata,
                    logx=False, logy=False,
                    xlabel_frac=0.8,
                    xlim=None,
                    debug=False,
                    pos="top right",
):
    """
    Plot linear regression of `ydata` on `xdata` using axes `ax`

    The plot is labelled with the regression slope (with uncertainty)
    and intercept, together with correlation coefficient, and p-value.

    If the axes are shown on a log scale, it may be desired to take
    log10 of data before regressing.  This is controlled by optional
    flag arguments `logx` and `logy`.

    Remaining arguments control aesthetics and layout of labeling
    """

    # Helper functions for optional log10 transformations of data
    def _x(x):
        "Forward data transformation for x axis"
        if logx:
            return np.log10(x)
        else:
            return x

    def _xx(x):
        "Backward data transformation for x axis"
        if logx:
            return 10**x
        else:
            return x

    def _y(y):
        "Forward data transformation for y axis"
        if logy:
            return np.log10(y)
        else:
            return y

    def _yy(y):
        "Backward data transformation for y axis"
        if logy:
            return 10**y
        else:
            return y


    # Linear regression on the residuals: y = a x + b
    # This does the Pearson r at the same time
    # a, b, r, p, da = linregress(_x(xdata), _y(ydata))
    a, b, r, p, da = robust_linregress(_x(xdata), _y(ydata))
    xstring = r"\log_{10} x" if logx else "x"
    ystring = r"\log_{10} y" if logy else "y"
    s = rf'Linear regression: ${ystring} = m \, {xstring} + c$'
    s += '\n' + f'$m = {a:.2f} \pm {da:.2f}$, $c = {b:.2f}$'
    s += '\n' + f'Correlation: $r = {r:.2f}$'
    if p > 0.001:
        s += f' $p = {p:.4f}$'
    else:
        pexp = np.floor(np.log10(p))
        s += rf' $p = {p/10**pexp:.1f} \times 10^{{{pexp:.0f}}}$'

    if debug:
        print(s)

    # xmin, xmax need to be specified in linear space
    if xlim is None:
        xlim = np.nanmin(xdata), np.nanmax(xdata)
    # Now transform to (maybe) log space
    _xgrid = np.linspace(_x(xlim[0]), _x(xlim[1]))
    # Center of x range for pivoting the slopes
    _xc = np.nanmean(_xgrid)
    # Where to place label, controlled by argument xlabel_frac
    _xlab = _xgrid[0] + xlabel_frac*(_xgrid[-1] - _xgrid[0])
    # We need to use the backward transform helper functions _xx and
    # _yy to get everything back to linear space for plotting (and
    # then let matplotlib decide whether to show results on log scale
    # or not)
    if debug:
        print("Arrow head position:", _xx(_xlab), _yy(a*_xlab + b))

    pos_options = "top left", "top right", "bottom left", "bottom right"
    assert pos in pos_options, "Invalid label position"
    vpos, hpos = pos.split()
    xytext = (
        0.02 if hpos == "left" else 0.98,
        1.05 if vpos == "top" else 0.02,
    )
    ha = hpos
    va = vpos

    ax.annotate(s=s, 
                xy=(_xx(_xlab), _yy(a*_xlab + b)),
                xytext=xytext, textcoords='axes fraction',
                arrowprops={'arrowstyle': '-|>', 'color': 'k'},
                fontsize='small', color='k',
                bbox={'fc': 'white', 'ec': 'none', 'alpha': 0.85, 'pad': 1},
                ha=ha, va=va)
    ax.fill_between(_xx(_xgrid),
                    _yy(a*_xc + (a+da)*(_xgrid - _xc) + b),
                    _yy(a*_xc + (a-da)*(_xgrid - _xc) + b),
                    alpha=0.4)
    ax.plot(_xx(_xgrid), _yy(a*_xc + a*(_xgrid - _xc) + b), lw=2, ls=':')
