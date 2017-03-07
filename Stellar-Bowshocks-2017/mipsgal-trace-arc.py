import glob
import os
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import astropy.coordinates as coord
from matplotlib import pyplot as plt
import seaborn as sns

sns.set_style('white')


SOURCE_DIR = 'OB/Kobulnicky2016'
source_table = Table.read(
    os.path.join(SOURCE_DIR, 'table1.dat'),
    format='ascii.cds',
    readme=os.path.join(SOURCE_DIR, 'ReadMe')
)

IMAGE_DIR = 'OB/MipsGal'

def skycoord_from_table_row(data):
    ra = f"{data['RAh']} {data['RAm']} {data['RAs']}"
    dec = f"{data['DE-']}{data['DEd']} {data['DEm']} {data['DEs']}"
    return coord.SkyCoord(f'{ra} {dec}', unit=(u.hourangle, u.deg))

for source_data in source_table[0:50]:
    # Coordinates of source central star
    c = skycoord_from_table_row(source_data)
    # Find all the images for this source
    provisional_list = glob.glob(
        f"{IMAGE_DIR}/*-{source_data['Name']}-*.fits")
    # Look for an image that is good
    good_image_list = []
    for image_name in provisional_list:
        hdu, = fits.open(image_name)
        looks_good = hdu.header['NAXIS1'] == hdu.header['NAXIS2']
        if looks_good:
            good_image_list.append(image_name)

    if good_image_list:
        # Use the first one in the list - because: why not? 
        hdu, = fits.open(good_image_list[0])
    else:
        # If there were no good images, then never mind
        continue

    # Create WCS object for this image
    w = WCS(hdu)
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
    rad_mask = (rpix > 0.1*R0) & (rpix < 3.0*R0)

    # Minimum and median brightness, which we might need later
    bright_min = np.nanmin(hdu.data)
    bright_median = np.nanmedian(hdu.data)

    # Next, a two step process to trace the arc.  First step is just
    # to find center and radius of curvature

    # Loop over a grid of angles between +/- 60 degrees
    ntheta = 25
    theta_grid, dtheta = np.linspace(-120.0, 120.0, ntheta, retstep=True)
    # Make everything be a longitude in range [-180:180]
    theta_grid = coord.Longitude(theta_grid, unit=u.degree, wrap_angle=180*u.degree)
    dtheta = coord.Longitude(dtheta, unit=u.degree, wrap_angle=180*u.degree)
    rpeak_grid = []
    rmean_grid = []
    bmax_grid = []
    bmean_grid = []
    for theta in theta_grid:
        # Select only pixels in the wedge within +/- dtheta/2 of this theta
        theta_mask = np.abs(theta_pix - theta) < 0.5*dtheta
        # Combine with the radius mask
        m = theta_mask & rad_mask

        if np.alltrue(~m):
            # If mask is empty, fill in this theta with NaNs
            rpeak = np.nan*u.deg
            rmean = np.nan*u.deg
            bright_max = np.nan
            bmean = np.nan
        else:            
            # Try a variety of methods for determining the arc radius
            # at this theta ...

            # Peak brightness
            ipeak = hdu.data[m].argmax()
            rpeak = rpix[m][ipeak]

            # Mean brightness-weighted radius. We divide weights by
            # radius to compensate for density of pixels. Also, we
            # select only points brighter than 0.5 times the peak
            # brightness in this wedge.  And all brightnesses are
            # relative to a background floor, which is either the
            # median over the whole image or, if the peak in the wedge
            # is lower than that, then the minimum over the image
            bright_max = hdu.data[m].max()
            bright_floor = bright_median if bright_max > bright_median else bright_min
            mb = (hdu.data - bright_floor) > 0.5*(bright_max - bright_floor)
            weights = (hdu.data[m & mb] - bright_floor)/rpix[m & mb]
            rmean = np.average(rpix[m & mb], weights=weights)*u.arcsec
            bmean = np.average(hdu.data[m & mb], weights=weights)

        # Fit Gaussian to profile - TODO?

        # Save all quantities into grid lists
        rmean_grid.append(rmean)
        rpeak_grid.append(rpeak)
        bmax_grid.append(bright_max)
        bmean_grid.append(bmean)

    # convert to single array of each quantity
    rmean_grid = coord.Angle(rmean_grid)
    rpeak_grid = coord.Angle(rpeak_grid)
    bmax_grid = np.array(bmax_grid)
    bmean_grid = np.array(bmean_grid)

    # Get the arc coordinates in RA, Dec
    #
    # Use the offset frame centered on the source and aligned with PA
    # axis.  The order of components is (lon, lat) where lon is
    # perpendicular and lat parallel to the bowshock axis
    rmean_coords = coord.SkyCoord(
        rmean_grid*np.sin(theta_grid),
        rmean_grid*np.cos(theta_grid),
        frame=offset_frame).transform_to('icrs')
    rpeak_coords = coord.SkyCoord(
        rpeak_grid*np.sin(theta_grid),
        rpeak_grid*np.cos(theta_grid),
        frame=offset_frame).transform_to('icrs')

    # Save a figure for each source
    fig = plt.figure(figsize=(10, 8))
    # Make a plot of the radii and brightnesses versus theta
    ax_r = plt.subplot(2, 2, 1)
    ax_i = plt.subplot(2, 2, 2, projection=w)
    ax_b = plt.subplot(2, 2, 3)
    ax_r.plot(theta_grid, rmean_grid, label='mean')
    ax_r.plot(theta_grid, rpeak_grid, label='peak')
    ax_r.axhline(R0.value)
    ax_r.legend()
    ax_r.set(ylim=[0.0, None], ylabel='Bow shock radius, arcsec')
    ax_b.plot(theta_grid, bmean_grid - bright_median, label='mean')
    ax_b.plot(theta_grid, bmax_grid - bright_median, label='peak')
    ax_b.legend()
    ax_b.set(ylim=[0.0, None], ylabel='Bow shock brightness',
             xlabel='Angle from nominal axis, degree'
    )

    # And also plot the image
    ax_i.imshow(hdu.data,
                vmin=bright_min, vmax=bmean_grid.max(), origin='lower')
    # Add a marker for the source
    wtran = ax_i.get_transform('world')
    ax_i.scatter(c.ra.deg, c.dec.deg, transform=wtran,
                 s=100, edgecolor='k', facecolor='orange')
    # And markers for the traced bow shock
    ax_i.scatter(rpeak_coords.ra.deg, rmean_coords.dec.deg, transform=wtran,
                 marker='+', c='r', s=30, alpha=0.5)

    fig.tight_layout()
    fig.savefig(image_name.replace('.fits', '-multiplot.png'))
