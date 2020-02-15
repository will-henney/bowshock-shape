import glob
import os
import sys
from pathlib import Path
from collections import OrderedDict
import numpy as np
from scipy.signal import medfilt, decimate
from astropy.table import Table, QTable
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import astropy.coordinates as coord
from astropy.modeling import models, fitting
from astropy.visualization.wcsaxes import SphericalCircle
from matplotlib import pyplot as plt
import seaborn as sns
import circle_fit_utils

sns.set_style('white')

DEBUG = False

SOURCE_DIR = Path('OB/Kobulnicky2016')
source_table = Table.read(
    str(SOURCE_DIR / 'table1.dat'),
    format='ascii.cds',
    readme=str(SOURCE_DIR / 'ReadMe')
)



##
## Deal with command-line arguments 
##

try:
    CIRCLE_THETA = abs(float(sys.argv[1]))*u.deg
    TRACE_METHOD = "peak" if float(sys.argv[1]) < 0 else "both"
except:
    print("Adopting default CIRCLE_THETA = 60")
    CIRCLE_THETA = 60.0*u.deg
    TRACE_METHOD = "both"

try:
    iseq_min, iseq_max = int(sys.argv[2]), int(sys.argv[3])
except:
    iseq_min, iseq_max = 0, 999_999_999

try:
    IMAGE_DIR = Path(sys.argv[4])
except:
    IMAGE_DIR = Path('OB/MipsGal')

try:
    imglob = sys.argv[5]
except:
    imglob = ""

# 10 Feb 2020 - now use a dedicated folder for output, which is named
# after the value of CIRCLE_THETA
OUTPUT_DIR = IMAGE_DIR / Path(f"circ_theta_{int(CIRCLE_THETA.value):03d}_{TRACE_METHOD}")

# Create the output folder if it does not already exist
OUTPUT_DIR.mkdir(exist_ok=True)

ENVIRONMENTS = {
    'I': 'Isolated',
    'H': 'H II region',
    'FH': 'Facing H II region',
    'FB': 'Facing bright-rimmed cloud',
}


# Some data in the the Kobulnicky2016 table is just wrong
OVERRIDE = {
    8: {'R0': 5.0},
    10: {'R0': 16.0},
    85: {'R0': 20.0},
    228: {'R0': 18.0},
    648: {'R0': 40},
    650: {'R0': 50},
}
STEP_BACK_FACTOR = 2.0

THMIN, THMAX = coord.Angle([-160.0*u.deg, 160.0*u.deg])

def description_from_table_row(data):
    desc = data['Name'] + '\n'
    if data['Alias']:
        desc += data['Alias'] + '\n'
    desc += f"R0 = {data['R0']:.1f} arcsec, PA = {data['PA']} deg" + '\n'
    csource = 'Multiple candidates' if data['Unc'] == 'C' else 'Single candidate'
    desc += f'{csource} for central source' + '\n'
    desc += f"Environment: {ENVIRONMENTS[data['Env']]}"
    return desc


def skycoord_from_table_row(data):
    ra = f"{data['RAh']} {data['RAm']} {data['RAs']}"
    dec = f"{data['DE-']}{data['DEd']} {data['DEm']} {data['DEs']}"
    return coord.SkyCoord(f'{ra} {dec}', unit=(u.hourangle, u.deg))


def coord_concat(array_tuple, **kwds):
    """Like numpy.concatenate but ensures result is of coordinate type"""
    return coord.SkyCoord(np.concatenate(array_tuple, **kwds))


def find_x_for_y_fraction(x, y, ymax, yfrac):
    """Find `x` where `y` is a fraction `yfrac` of maximum `ymax`

    Returns tuple of left and right x values.  Assumes unimodal peak
    in y(x) and that the peak is near to x = 0.  If y never falls
    below the fraction on the left (right), then returns the minimum
    (maximum) value of x on that side.
    """
    try:
        xleft = np.max(x[(y < yfrac*ymax) & (x < 0.0)])
    except:
        xleft = x.min()
    try:
        xright = np.min(x[(y < yfrac*ymax) & (x > 0.0)])
    except:
        xright = x.max()
    return xleft, xright


fitter = fitting.LevMarLSQFitter()

#
# Sherpa/Saba is suspended 11 Feb 2020 - crashes with composite models
#
# This one is for getting uncertainties in the fitted shell thickness
# It require sherpa and saba packages to be installed (I had to
# install them from source)
# sfitter = fitting.SherpaFitter(
#     statistic='chi2',
#     optimizer='levmar',
#     estmethod='confidence',
# )

for source_data in source_table:

    if not iseq_min <= source_data['Seq'] <= iseq_max:
        # Skip this source if not in the requested range
        continue

    print(source_data)

    # Override data from table where necessary
    if source_data['Seq'] in OVERRIDE:
        for k, v in OVERRIDE[source_data['Seq']].items():
            source_data[k] = v

    # Coordinates of source central star
    c = skycoord_from_table_row(source_data)
    # Find all the images for this source
    provisional_list = IMAGE_DIR.glob(
        # f"*-{source_data['Name']}-*.fits"
        f"{source_data['Seq']:04d}*{imglob}.fits"
    )
    # Look for an image that is good
    good_image_list = []
    for image_path in provisional_list:
        hdu, = fits.open(image_path)
        # looks_good = hdu.header['NAXIS1'] == hdu.header['NAXIS2']
        looks_good = (hdu.header['NAXIS1'] > 100) &  (hdu.header['NAXIS2'] > 100)
        if looks_good:
            good_image_list.append(image_path)

    if good_image_list:
        # Use the first one in the list - because: why not?
        image_path = good_image_list[0]
        hdu, = fits.open(image_path)
    else:
        # If there were no good images, then never mind
        print("ABORT - no good images")
        continue

    image_stem = image_path.stem

    # Create WCS object for this image
    w = WCS(hdu)
    # Save pixel scale for later use
    pix_size = coord.Angle(abs(w.wcs.cdelt[0]), unit=u.deg)
    # Find celestial coordinates for each pixel
    ny, nx = hdu.data.shape
    xpix = np.arange(nx)[None, :]
    ypix = np.arange(ny)[:, None]
    cpix = coord.SkyCoord.from_pixel(xpix, ypix, w)
    # Now find radius and position angle from source
    rpix = c.separation(cpix).to(u.arcsec)
    pa_pix = c.position_angle(cpix).to(u.degree)
    # Nominal PA of bowshock axis from table
    pa0 = coord.Angle(source_data['PA'], unit=u.degree)
    # theta is angle from nominal axis, set to range [-180:180]
    theta_pix = coord.Longitude(pa_pix - pa0, wrap_angle=180*u.degree)

    # Also create an offset frame in case we need it later, in which
    # measurements are with respect to the central source coordinate,
    # and rotated by pa0
    offset_frame = c.skyoffset_frame(rotation=pa0)
    # Coordinates of each pixel in the offset frame - this has
    # components: offset_pix.lat (along pa0) and offset_pix.lon
    # (perpendicular to pa0)
    offset_pix = cpix.transform_to(offset_frame)

    # Nominal arc radius from source table
    R0 = source_data['R0']*u.arcsec

    # Only look in a restricted range of radius around R0
    rad_mask = (rpix > 0.5*R0) & (rpix < 3.0*R0)

    # Minimum and median brightness, which we might need later
    bright_min = np.nanmin(hdu.data)
    bright_median = np.nanmedian(hdu.data)
    # Calculate a robust stddev as 1.4826 times MAD
    bright_mad = np.nanmedian(np.abs(hdu.data - bright_median))
    bright_sigma = 1.4826*bright_mad

    # Next, we trace the arc

    # 08 Mar 2017 - Try a different tack - take radii from a center
    # that is "stepped back" by STEP_BACK_FACTOR times R0 away from
    # the source
    c_sb = coord.SkyCoord(0.0*u.deg, -STEP_BACK_FACTOR*R0,
                          frame=offset_frame).transform_to('icrs')
    # Repeat all the above to find radius, angle from this new point
    # Now find radius and position angle from source
    r_sb_pix = c_sb.separation(cpix).to(u.arcsec)
    pa_sb_pix = c_sb.position_angle(cpix).to(u.degree)
    # theta is angle from nominal axis, set to range [-180:180]
    th_sb_pix = coord.Longitude(pa_sb_pix - pa0, wrap_angle=180*u.degree)
    # And a frame relative to the "step back" center too
    sb_offset_frame = c_sb.skyoffset_frame(rotation=pa0)


    # Loop over a grid of angles between +/- 80 degrees
    ntheta = 68
    theta_grid, dtheta = np.linspace(-80.0, 80.0, ntheta, retstep=True)
    # Make everything be a longitude in range [-180:180]
    th_sb_grid = coord.Longitude(theta_grid, unit=u.degree, wrap_angle=180*u.degree)
    dtheta = coord.Longitude(dtheta, unit=u.degree, wrap_angle=180*u.degree)
    r_sb_peak_grid = []
    r_sb_mean_grid = []
    bmax_grid = []
    bmean_grid = []
    for th_sb in th_sb_grid:
        # Select only pixels in the wedge within +/- dtheta/2 of this theta
        theta_mask = np.abs(th_sb_pix - th_sb) < 0.5*dtheta
        # Combine with the radius mask
        m = theta_mask & rad_mask

        if np.alltrue(~m):
            # If mask is empty, fill in this theta with NaNs
            r_sb_peak = np.nan*u.deg
            r_sb_mean = np.nan*u.deg
            bright_max = np.nan
            bmean = np.nan
            print("Empty mask for theta =", th_sb)
        else:            
            # Try a variety of methods for determining the arc radius
            # at this theta ...

            # Peak brightness
            ipeak = hdu.data[m].argmax()
            r_sb_peak = r_sb_pix[m][ipeak]

            # Mean brightness-weighted radius. We divide weights by
            # radius to compensate for density of pixels. Also, we
            # select only points brighter than 0.5 times the peak
            # brightness in this wedge.  And all brightnesses are
            # relative to a background floor, which is either the
            # median over the whole image or, if the peak in the wedge
            # is lower than that, then the minimum over the image
            bright_max = np.nanmax(hdu.data[m])
            bright_floor = bright_median if bright_max > bright_median else bright_min
            mb = (hdu.data - bright_floor) > 0.5*(bright_max - bright_floor)
            weights = (hdu.data[m & mb] - bright_floor)/r_sb_pix[m & mb]
            try: 
                r_sb_mean = np.average(r_sb_pix[m & mb], weights=weights)
                bmean = np.average(hdu.data[m & mb], weights=weights)
            except ZeroDivisionError:
                r_sb_mean = np.nan*u.deg
                bmean = np.nan
            if DEBUG:
                print("Bright max, mean = ", bright_max, bmean, "for theta =", th_sb)

        # Fit Gaussian to profile - TODO?

        # Save all quantities into grid lists
        r_sb_mean_grid.append(r_sb_mean)
        r_sb_peak_grid.append(r_sb_peak)
        bmax_grid.append(bright_max)
        bmean_grid.append(bmean)

    # convert to single array of each quantity
    r_sb_mean_grid = coord.Angle(r_sb_mean_grid)
    r_sb_peak_grid = coord.Angle(r_sb_peak_grid)
    bmax_grid = np.array(bmax_grid)
    bmean_grid = np.array(bmean_grid)


    # Now switch back to the frame centered on the source

    # Get the arc coordinates in RA, Dec
    #
    # Use the offset frame centered on the source and aligned with PA
    # axis.  The order of components is (lon, lat) where lon is
    # perpendicular and lat parallel to the bowshock axis
    rmean_coords = coord.SkyCoord(
        r_sb_mean_grid*np.sin(th_sb_grid),
        r_sb_mean_grid*np.cos(th_sb_grid),
        frame=sb_offset_frame).transform_to('icrs')
    rpeak_coords = coord.SkyCoord(
        r_sb_peak_grid*np.sin(th_sb_grid),
        r_sb_peak_grid*np.cos(th_sb_grid),
        frame=sb_offset_frame).transform_to('icrs')

    # Switch back to frame centered on source
    rmean_grid = c.separation(rmean_coords).to(u.arcsec)
    theta_mean_grid = coord.Longitude(
        c.position_angle(rmean_coords).to(u.degree) - pa0,
        wrap_angle=180*u.degree)
    rpeak_grid = c.separation(rpeak_coords).to(u.arcsec)
    theta_peak_grid = coord.Longitude(
    c.position_angle(rpeak_coords).to(u.degree) - pa0,
        wrap_angle=180*u.degree)

    # Fit circle to peak points within CIRCLE_THETA of axis
    cmask_peak = np.abs(theta_peak_grid) <= CIRCLE_THETA
    cmask_mean = np.abs(theta_mean_grid) <= CIRCLE_THETA
    # Use the mean and peak points
    try:
        if TRACE_METHOD == "both":
            points2fit = coord_concat((rpeak_coords[cmask_peak],
                                       rmean_coords[cmask_mean]))
        elif TRACE_METHOD == "peak":
            points2fit = rpeak_coords[cmask_peak]

        # Winnow out the points with radii too far from median
        r2fit = c.separation(points2fit).to(u.arcsec)
        rmed = np.median(r2fit)
        print('Median radius for circle fit:', rmed)
        mfit = (r2fit >= 0.5*rmed) & (r2fit <= 2.0*rmed)
        n_drop = (~mfit).sum()
        if n_drop > 0:
            print(n_drop, 'points dropped for circle fit')
        points2fit = points2fit[mfit]
        # Initial guess for center would make Rc/R0 = 2
        center0 = coord.SkyCoord(0.0*u.deg, -R0,
                                 frame=offset_frame).transform_to('icrs')
        Rc, center = circle_fit_utils.fit_circle(points2fit, center0)
        Rc = Rc.to(u.arcsec)
    except:
        print('ABORT - Problem with points2fit or circle_fit_utils.fit_circle')
        continue


    if Rc > 100*R0:
        print('ABORT due to ridiculous radius of curvature: Rc =', Rc)
        continue

    # Find standard deviation of points from circle
    Rc_sigma = np.std(
        circle_fit_utils.deviation_from_circle(points2fit, center)
    ).to(u.arcsec)

    # Find PA of circle fit
    # First assume case where center of curvature is "behind" the source
    pa_circ = center.position_angle(c).to(u.deg)
    # Find difference between fitted and nominal position angle
    delta_pa = coord.Longitude(pa_circ - pa0, wrap_angle=180*u.deg)
    if np.abs(delta_pa) > 90*u.deg:
        # Check for Case where center of curvature is "in front of" the source
        pa_circ = c.position_angle(center).to(u.deg)
        delta_pa = coord.Longitude(pa_circ - pa0, wrap_angle=180*u.deg)

    # Find our estimate of R0
    #
    # Make some masks selecting points within 10 deg of pa_circ
    # m0_peak = np.abs(theta_peak_grid - delta_pa) <= 10.0*u.deg
    # m0_mean = np.abs(theta_mean_grid - delta_pa) <= 10.0*u.deg
    # 05 Feb 2020 - change to using nominal PA, instead of PA_circ to define axis
    m0_peak = np.abs(theta_peak_grid) <= 10.0*u.deg
    m0_mean = np.abs(theta_mean_grid) <= 10.0*u.deg

    # Then concatenate all the R values that meet this condition
    if TRACE_METHOD == "both":
        R0_grid = coord.Angle(
            np.concatenate((rpeak_grid.value[m0_peak],
                            rmean_grid.value[m0_mean])),
            unit=u.arcsec)
    elif TRACE_METHOD == "peak":
        R0_grid = coord.Angle(rpeak_grid.value[m0_peak], unit=u.arcsec)

    R0_fit, R0_sigma = R0_grid.mean(), R0_grid.std()
    fit_msg = f'Fitted R0 = {R0_fit.arcsec:.1f} +/- {R0_sigma.arcsec:.1f} arcsec' + '\n'

    # Make an offset frame centered on the center of curvature
    circ_offset_frame = center.skyoffset_frame(rotation=pa_circ)
    # Find R(theta) for the fitted circle
    thdash = np.linspace(-180.0, 180.0, 501)*u.deg
    circ_points = coord.SkyCoord(
        Rc*np.sin(thdash), Rc*np.cos(thdash),
        frame=circ_offset_frame).transform_to('icrs')
    circ_theta = coord.Longitude(
        c.position_angle(circ_points).to(u.deg) - pa0,
        wrap_angle=180*u.degree)
    circ_radius = c.separation(circ_points).to(u.arcsec)
    # Eliminate points that are not within +/- 100 deg of nominal axis
    mcirc = (circ_theta >= -100.0*u.deg) & (circ_theta <= 100.0*u.deg)
    circ_radius[~mcirc] *= np.nan
    fit_msg += f'PA_circ = {pa_circ.deg:.1f}, delta PA = {delta_pa.deg:.1f}' + '\n'
    fit_msg += f'Rc = {Rc.arcsec:.1f} +/- {Rc_sigma.arcsec:.1f} arcsec' + '\n'

    # Find R90
    #
    # Make some masks selecting points within 10 deg of +90 and -90
    m90p_peak = np.abs(theta_peak_grid - 90.0*u.deg) <= 10.0*u.deg
    m90n_peak = np.abs(theta_peak_grid + 90.0*u.deg) <= 10.0*u.deg
    m90p_mean = np.abs(theta_mean_grid - 90.0*u.deg) <= 10.0*u.deg
    m90n_mean = np.abs(theta_mean_grid + 90.0*u.deg) <= 10.0*u.deg
    # Then concatenate all the R values for the two cases
    if TRACE_METHOD == "both":
        R90p_grid = coord.Angle(
            np.concatenate((rpeak_grid.value[m90p_peak],
                            rmean_grid.value[m90p_mean])),
            unit=u.arcsec)
        R90n_grid = coord.Angle(
            np.concatenate((rpeak_grid.value[m90n_peak],
                            rmean_grid.value[m90n_mean])),
            unit=u.arcsec)
    elif TRACE_METHOD == "peak":
        R90p_grid = coord.Angle(rpeak_grid.value[m90p_peak], unit=u.arcsec)
        R90n_grid = coord.Angle(rpeak_grid.value[m90n_peak], unit=u.arcsec)

    # And calculate mean and standard deviation
    R90p, R90p_sigma = R90p_grid.mean(), R90p_grid.std()
    R90n, R90n_sigma = R90n_grid.mean(), R90n_grid.std()
    fit_msg += f'R90+ = {R90p.arcsec:.1f} +/- {R90p_sigma.arcsec:.1f} arcsec' + '\n'
    fit_msg += f'R90- = {R90n.arcsec:.1f} +/- {R90n_sigma.arcsec:.1f} arcsec' + '\n'

    # Calculate angular fwhm width of brightness distribution
    x, y = theta_peak_grid.value, bmax_grid - bright_median
    try:
        b_max = np.nanmax(bmax_grid[m0_peak]) - bright_median
    except ValueError:
        print('ABORT - No valid pixels in m0_peak, cannot continue')
        continue

    if not np.isfinite(b_max):
        print("Axis peak brightnesses:", bmax_grid[m0_peak])
        print("Setting to 10 x median:", bright_median)
        b_max = 20*bright_median

    # Use 3 different brightness levels: [0.841, 0.5, 0.210], which
    # for a Gaussian correspond to [0.5, 1, 1.5] times the FWHM
    th841_left, th841_right = find_x_for_y_fraction(x, y, b_max, 0.841)
    th500_left, th500_right = find_x_for_y_fraction(x, y, b_max, 0.5)
    th210_left, th210_right = find_x_for_y_fraction(x, y, b_max, 0.210)

    # Another attempt at angular half width using Gaussian fit
    g_init = (models.Gaussian1D(amplitude=b_max, mean=0.0, stddev=60.0)
              + models.Const1D(amplitude=0.05*b_max))
    g_fit = fitter(g_init, x, y)

    # Finally fit the width on the axis

    # THE RADIAL BRIGHTNESS PROFILE
    #
    # This is to try and measure h, the bow thickness
    #
    # Use the offset_pix frame, centered on the star, with
    # offset_pix.lat along the axis and offset_pix.lon across the axis.

    # Restrict to a strip of width 0.3 times fitted radius
    strip_breadth = 0.3*R0_fit
    m = np.abs(offset_pix.lon) < 0.5*strip_breadth
    # Estimate number of image pixels across the strip
    n_across_strip = int(round(strip_breadth.arcsec / pix_size.arcsec))

    # Also restrict to points not too far from star
    m = m & (rpix < 4.0*R0_fit) 
    x, y, across = offset_pix.lat[m].arcsec, hdu.data[m], offset_pix.lon[m].arcsec

    gg_init = (
        models.Gaussian1D(amplitude=np.nanmax(y),
                          mean=R0_fit.arcsec,
                          stddev=0.5*R0_fit.arcsec)
        + models.Gaussian1D(amplitude=0.1*np.nanmax(y),
                            mean=0.0,
                            stddev=6.0)
        + models.Polynomial1D(degree=2, c0=bright_median)
        )
    # Obligatory to provide errors with SherpaFitter, at least for "chi2" statistic
    # gg_fit = sfitter(gg_init, x, y, err=bright_median*np.ones_like(y))
    try:
        gg_fit = fitter(gg_init, x, y)
    except:
        print("ABORT: Failed to fit radial brightness profile")
        continue

    # 11 Feb 2020 - shelve Sherpa for now
    # # Find the 1-sigma confidence intervals
    # gg_errors = sfitter.est_errors(sigma=1.0)

    # Save derived values from fitting the brightness profile
    R0_from_bright_fit = gg_fit[0].mean.value
    # Multiply sigma by 2 sqrt(2 ln(2)) to get FWHM
    H_from_bright_fit = 2.35*gg_fit[0].stddev.value
    Peak_from_bright_fit = gg_fit[0].amplitude.value
    # Note we need to unpack a length-1 array here
    BG_from_bright_fit, = gg_fit(R0_from_bright_fit) - Peak_from_bright_fit
    Contrast_fom_bright_fit = Peak_from_bright_fit / BG_from_bright_fit

    # Alternative thickness measurement II: direct FWHM
    # This is the same as we did for the angle breadth

    # First subtract polynomial background
    yy = y - gg_fit[-1](x)
    # Center x axis on peak position
    xx = x - R0_fit.arcsec
    # Sort the x, y points by x
    pts_order = np.argsort(xx)
    xx = xx[pts_order]
    yy = yy[pts_order]
    # Apply median filter by typical number of points in the
    # transverse direction.
    n_across_strip = 2*(n_across_strip // 2) + 1  # ensure odd number 
    yy = medfilt(yy, n_across_strip)
    xx = medfilt(xx, n_across_strip)
    # Then winnow down the number of points by a similar factor
    nwinnow = max(1, n_across_strip // 2)
    yy = yy[::nwinnow]
    xx = xx[::nwinnow]
    yymax = np.nanmax(yy[np.abs(xx <= R0_fit.arcsec)])
    # Find Half-max points
    r_in, r_out = find_x_for_y_fraction(xx, yy, yymax, 0.5)
    H_direct_fwhm = r_out - r_in

    # Save the fit data for each source
    table_file_name = f"{image_stem}-arcfit.tab"
    save_vars = [
        ['Seq', source_data['Seq']], 
        ['R0_fit', R0_fit.arcsec], 
        ['R0_sigma', R0_sigma.arcsec],
        ['pa_circ', pa_circ.deg],
        ['delta_pa', delta_pa.deg],
        ['Rc', Rc.arcsec],
        ['Rc_sigma', Rc_sigma.arcsec],
        ['R90p', R90p.arcsec],
        ['R90p_sigma', R90p_sigma.arcsec],
        ['R90n', R90n.arcsec],
        ['R90n_sigma', R90n_sigma.arcsec],
        ['th_p', th500_right],
        ['th_m', th500_left],
        ['th_p_841', th841_right],
        ['th_m_841', th841_left],
        ['th_p_210', th210_right],
        ['th_m_210', th210_left],
        ['th_g', g_fit.mean_0.value],
        ['dth_g', 2.35*g_fit.stddev_0.value],
        ['H_g', H_from_bright_fit],
        ['H_d', H_direct_fwhm],
        ['H_d_out', r_out],
        ['H_d_in', r_in],
        ['R0_g', R0_from_bright_fit],
        ['Peak24', Peak_from_bright_fit],
        ['Contrast', Contrast_fom_bright_fit],
    ]
    colnames, colvals = zip(*save_vars)
    Table(rows=[list(colvals)],
          names=list(colnames)).write(
              OUTPUT_DIR / table_file_name,
              format='ascii.tab',
              overwrite=True)

    # Save a figure for each source
    fig = plt.figure(figsize=(12, 8))

    # Make a plot of the radii and brightnesses versus theta
    ax_r = fig.add_axes((0.08, 0.55, 0.35, 0.4))
    ax_b = fig.add_axes((0.08, 0.08, 0.35, 0.4))
    ax_i = fig.add_axes((0.5, 0.1, 0.45, 0.45), projection=w)
    ax_r.plot(theta_mean_grid, rmean_grid, 'o', c='c', label='mean')
    ax_r.plot(theta_peak_grid, rpeak_grid, 'o', c='r', label='peak')
    ax_r.plot(circ_theta, circ_radius,
              '--', c='m', label='circle fit')
    ax_r.axhline(R0.value)
    ax_r.axvspan(-90.0, 90.0, facecolor='k', alpha=0.05)
    ax_r.axvspan(-CIRCLE_THETA.value, CIRCLE_THETA.value, facecolor='k', alpha=0.05)
    ax_r.axvline(0.0, c='k', ls='--')
    ax_r.legend()
    ax_r.set(xlim=[THMIN.deg, THMAX.deg],
             ylim=[0.0, None], ylabel='Bow shock radius, arcsec')
    ax_b.plot(theta_mean_grid, bmean_grid - bright_median, 'o', c='c', label='mean')
    ax_b.plot(theta_peak_grid, bmax_grid - bright_median, 'o', c='r', label='peak')
    ax_b.plot([th500_left, th500_right], [0.5*b_max]*2,
              "-", c="r", alpha=0.2, lw=5, label="_nolabel_")
    ax_b.plot([th841_left, th841_right], [0.841*b_max]*2,
              "-", c="r", alpha=0.2, lw=3, label="_nolabel_")
    ax_b.plot([th210_left, th210_right], [0.210*b_max]*2,
              "-", c="r", alpha=0.2, lw=3, label="_nolabel_")
    thgrid = np.linspace(THMIN.deg, THMAX.deg, 101)
    ax_b.plot(thgrid, g_fit(thgrid), '--', c="r", alpha=0.3, lw=3, label="_nolabel_")
    ax_b.axvspan(-90.0, 90.0, facecolor='k', alpha=0.05)
    ax_b.axvspan(-CIRCLE_THETA.value, CIRCLE_THETA.value, facecolor='k', alpha=0.05)
    ax_b.axvline(0.0, c='k', ls='--')
    ax_b.legend()
    ax_b.set(xlim=[THMIN.deg, THMAX.deg], xlabel='Angle from nominal axis, degree',
             ylim=[0.0, 1.4*b_max],
             ylabel='Bow shock brightness',
    )

    # And also plot the image
    ax_i.imshow(hdu.data,
                vmin=bright_min, vmax=bmean_grid[m0_peak].max(), origin='lower')

    # And contours
    if bmax_grid[m0_peak].max() > bright_median:
        clevels = np.linspace(bright_median, bmax_grid[m0_peak].max(), 10)
    else:
        clevels = np.linspace(bright_median, hdu.data.max(), 10)
    ax_i.contour(hdu.data, levels=clevels, alpha=0.5)

    wtran = ax_i.get_transform('world')

    # Add markers for the traced bow shock
    ax_i.scatter(rmean_coords.ra.deg, rmean_coords.dec.deg, transform=wtran,
                 marker='.', c='c', s=30, alpha=0.5)
    ax_i.scatter(rpeak_coords.ra.deg, rpeak_coords.dec.deg, transform=wtran,
                 marker='.', c='r', s=30, alpha=0.5)

    # Add a line for the PA orientation
    PA_coords = coord.SkyCoord(
        [0.0*u.deg, 0.0*u.deg], [-2*R0, 2*R0],
        frame=offset_frame).transform_to('icrs')
    ax_i.plot(PA_coords.ra.deg, PA_coords.dec.deg,
              transform=wtran, c='orange', lw=2, alpha=0.8)

    # And plot the fitted circle
    circ = SphericalCircle((center.ra, center.dec), Rc,
                           edgecolor='m', lw=2, alpha=0.5, facecolor='none',
                           transform=wtran)
    ax_i.add_patch(circ)
    # And a line for the fitted PA axis
    PA_fit_coords = coord.SkyCoord(
        [0.0*u.deg, 0.0*u.deg], [-1.2*Rc, 1.2*Rc],
        frame=circ_offset_frame).transform_to('icrs')
    ax_i.plot(PA_fit_coords.ra.deg, PA_fit_coords.dec.deg,
              transform=wtran, c='m', lw=1.5, alpha=0.6)
    # And the center of curvature
    ax_i.scatter(center.ra.deg, center.dec.deg, transform=wtran,
                 marker='o', s=30, edgecolor='k', facecolor='m')
    # Add a marker for the source
    ax_i.scatter(c.ra.deg, c.dec.deg, transform=wtran,
                 s=150, marker='*', edgecolor='k', facecolor='orange')


    # Add coordinate grids
    ax_i.coords.grid(color='m', linestyle='solid', alpha=0.2)
    ax_i.coords['ra'].set_axislabel('Right Ascension')
    ax_i.coords['dec'].set_axislabel('Declination')
    overlay = ax_i.get_coords_overlay('galactic')
    overlay.grid(color='c', linestyle='solid', alpha=0.2)
    overlay['l'].set_axislabel('Galactic Longitude')
    overlay['b'].set_axislabel('Galactic Latitude')

    # Add title
    ax_i.text(0.5, 1.7, description_from_table_row(source_data),
              transform=ax_i.transAxes, ha='center', va='bottom'
    )
    ax_i.text(0.5, 1.6, fit_msg,
              transform=ax_i.transAxes, ha='center', va='top'
    )

    fig.savefig(OUTPUT_DIR /  f"{image_stem}-multiplot.pdf")
    # Important to close figure explicitly so as not to leak resources
    plt.close(fig)



    # And another fig for the radial brightness profile
    fig, ax = plt.subplots()
    ax.scatter(x, y, c=across, marker='.', alpha=0.2)
    xgrid = np.linspace(x.min(), x.max(), 101)
    ax.plot(xgrid, gg_fit(xgrid), '-', color="k")
    ax.plot(xgrid, gg_fit[0](xgrid), '--', lw=0.5, color="k")
    # And plot the median-filtered profile used to find direct FWHM
    ax.plot(xx + R0_fit.arcsec, yy, '-', color="g")
    ax.plot(np.array([r_in, r_out]) + R0_fit.arcsec, [0.5*yymax]*2,
            "-", color="g", lw=3, alpha=0.3)
    ax.axvline(0.0, c="k", ls="--")
    ax.axvline(R0_fit.arcsec, c="k", ls=":")
    ax.set(
        xlabel = "Offset along axis, arcsec",
        ylabel = r"Brightness in strip of width $0.3 \times R_0$",
    )
    fig.savefig(OUTPUT_DIR / f'{image_stem}-rad-bright.pdf')
    plt.close(fig)
