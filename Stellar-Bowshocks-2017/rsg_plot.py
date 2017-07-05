import json
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization.lupton_rgb import LinearMapping
from reproject import reproject_interp

def plot_rsg_image_fit(source='AlphaOri', batch=1, 
                       pmin=1, pmax=99,
                       pmin160=1, pmax160=99,
                       extra='', extra160='',
                       window_factor=1.5,
                       show160=False,
                       label=r'$\alpha$ Ori'):
    jsonfile = f'RSG/{source.lower()}-arcdata.json'
    with open(jsonfile) as f:
        arcdata = json.load(f)

    fits_pattern = (f'RSG/MESS_PPHOT/Batch{batch:02d}/'
                    + f'{source.upper()}*{extra}_70_*pixfrac*.fits')
    fitsfile = glob.glob(fits_pattern)[0]
    hdulist = fits.open(fitsfile)
    hdu = hdulist['image']
    hdu160 = fits.open(fitsfile.replace(f'{extra}_70',
                                        f'{extra160}_160'))['image']

    # Fix the unit keywords
    for k in 'CUNIT1', 'CUNIT2':
        hdu.header[k] = 'deg'
        hdu160.header[k] = 'deg'

    # Reproject the 160 image onto the same grid as 70
    hdu160_remap, mask = reproject_interp(hdu160, hdu.header)

    w = WCS(hdu.header)
    xpix_arcsec, ypix_arcsec = 3600*w.wcs.cdelt

    # Find plot window limits
    xrpix, yrpix = w.wcs.crpix
    window_size = window_factor*arcdata['outer']['Rc']
    x1, x2 = xrpix + window_size/xpix_arcsec, xrpix - window_size/xpix_arcsec
    y1, y2 = yrpix - window_size/ypix_arcsec, yrpix + window_size/ypix_arcsec

    # Two axes side by side
    fig, (axim, ax) = plt.subplots(1, 2, figsize=(10, 5),
                                   sharex=True, sharey=True,
                                   subplot_kw=dict(projection=w))

    # Brightness limits for images
    vmin, vmax = np.percentile(hdu.data, [pmin, pmax])
    vmin160, vmax160 = np.percentile(hdu160.data, [pmin160, pmax160])

    # Make an RGB image from the 80 and 160 bands
    red = (hdu160_remap - vmin160)/(vmax160 - vmin160)
    blue = (hdu.data - vmin)/(vmax - vmin)
    green = (red + blue)/2.0
    rgbmap = LinearMapping(minimum=0.0, maximum=1.0)
    rgbim = rgbmap.make_rgb_image(red, green, blue)

    axim.imshow(rgbim)
    # axim.imshow(hdu.data, origin='lower',
    #             vmin=vmin, vmax=vmax, cmap='gray', alpha=1.0)
    if show160:
        ax.imshow(hdu160_remap, origin='lower',
                  vmin=vmin160, vmax=vmax160, cmap='gray_r', alpha=0.6)
    else:
        ax.imshow(hdu.data, origin='lower',
                  vmin=vmin, vmax=vmax, cmap='gray_r', alpha=0.6)
    bmax = hdu.data.max()
    axim.contour(hdu.data, levels=[0.05*bmax, 0.2*bmax], colors='r')
    ax.contour(hdu.data, levels=[0.05*bmax, 0.2*bmax], colors='r')

    # Star position
    [x0], [y0] = w.all_world2pix([arcdata['star']['RA_dg']],
                                 [arcdata['star']['Dec_dg']], 0)

    # Plot circle fit
    xc0 = x0 + arcdata['outer']['xc']/xpix_arcsec
    yc0 = y0 + arcdata['outer']['yc']/ypix_arcsec
    ax.scatter(xc0, yc0, marker='o', s=30, edgecolor='k', facecolor='m',
               zorder=5)
    circ = matplotlib.patches.Circle((xc0, yc0),
                                     radius=arcdata['outer']['Rc']/ypix_arcsec,
                                     edgecolor='m', lw=2, alpha=0.5, facecolor='none',
    )
    ax.add_patch(circ)

    # Plot fitted PA axis
    cvec = np.array([(xc0 - x0), (yc0 - y0)])
    cvec /= np.sqrt((xc0 - x0)**2 + (yc0 - y0)**2)
    cvec *= 1.2*arcdata['outer']['Rc']/ypix_arcsec
    ax.plot([xc0 + cvec[0], xc0 - cvec[0]],
            [yc0 + cvec[1], yc0 - cvec[1]],
            c='m', lw=1.5, alpha=0.6)

    # Plot central star
    ax.scatter(x0, y0, s=150, marker='*', edgecolor='k', facecolor='orange', zorder=6)

    # Plot arc points
    xp = x0 + arcdata['outer']['x']/xpix_arcsec
    yp = y0 + arcdata['outer']['y']/ypix_arcsec
    ax.scatter(xp, yp, marker='o', s=10, c='r', alpha=0.3, zorder=10)

    # Add label with source name and R0
    R0_arcmin = arcdata['outer']['R0']/60
    s = rf"{label}    $R_0 = {R0_arcmin:.2f}'$"
    whitebox = {'fc': 'white', 'alpha': 0.5}
    ax.text(0.05, 0.05, s,
            ha='left', va='baseline', bbox=whitebox,
            transform=ax.transAxes)

    ax.set(xlim=[x1, x2], ylim=[y1, y2])
    axim.set(xlim=[x1, x2], ylim=[y1, y2])
    axim.coords['ra'].set_axislabel('Right Ascension')
    axim.coords['dec'].set_axislabel('Declination')
    ax.coords['dec'].set_ticklabel_visible(False)
    figfile = f'RSG/{source.lower()}-imageplot.pdf'
    fig.savefig(figfile)
    plt.close(fig)
    return figfile
