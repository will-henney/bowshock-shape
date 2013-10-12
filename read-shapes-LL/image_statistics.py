"""
Library of functions for calculating statistics

Some originally came from hires-extract/measure-noise.py
"""

import numpy as np
import scipy.special as sp


def trimean_and_iqr(data):
    # Find the 1st, 2nd, 3rd quartiles
    q1 = np.percentile(data, 25.0)
    q2 = np.percentile(data, 50.0)
    q3 = np.percentile(data, 75.0)
    # The trimean is the weighted average of the median and the
    # two quartiles (http://en.wikipedia.org/wiki/Trimean)
    trimean = 0.25*(q1 + 2*q2 + q3)
    # The interquartile range is the difference between quartiles
    iqr = q3 - q1
    return trimean, iqr


def robust_statistics(x, y, xedges):
    """Calculate robust estimates of location and scale of a distribution

    Returns (loc, scale) of y, binned by x according to xedges

    Returns vectors of length len(xedges)-1 

    loc is the "average" value, estimated as the trimean

    scale is the width of the distribution, estimated from the
    interquartile range (IQR) and rescaled to be equal to the standard
    deviation, sigma, for a Gaussian distribution

    The IQR and trimean are much less affected by outliers and
    power-law wings than are the regular mean and sigma

    """
    # Interquartile range in units of std for Gaussian profile
    iqr_over_sigma = 2.0*np.sqrt(2.0)*sp.erfinv(0.5)
    trimeans, iqrs = [], []
    # Loop over bins in x
    for x1, x2 in zip(xedges[:-1], xedges[1:]):
        # Construct a mask that selects this bin
        m = (x >= x1) & (x < x2)
        if m.sum() > 1:
            trimean, iqr = trimean_and_iqr(y[m])
        else:
            trimean, iqr = np.nan, np.nan
    
        trimeans.append(trimean)
        iqrs.append(iqr)
    # Convert lists to numpy arrays before returning
    loc = np.array(trimeans)
    # Put scale in units of Gaussian standard deviation
    scale = np.array(iqrs)/iqr_over_sigma
    return loc, scale 
