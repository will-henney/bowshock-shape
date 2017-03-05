import os
import requests
import xmltodict
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import astropy.coordinates as coord

SST_URL = 'http://irsa.ipac.caltech.edu/cgi-bin/Atlas/nph-atlas'
mipsgal_params = {
    'mission': 'MIPSGAL',
    'mode': 'PI',
    'regSize': '0.01',
    'covers': 'on',
}
IMG_URL_ROOT = 'https://irsa.ipac.caltech.edu:443/data/SPITZER/MIPSGAL'

SOURCE_DIR = 'OB/Kobulnicky2016'
source_table = Table.read(
    os.path.join(SOURCE_DIR, 'table1.dat'),
    format='ascii.cds',
    readme=os.path.join(SOURCE_DIR, 'ReadMe')
)

OUTPUT_IMAGE_DIR = 'OB/MipsGal'
IMAGE_SIZE_DEGREES = 4.0/60.0           

def skycoord_from_table_row(data):
    ra = f"{data['RAh']} {data['RAm']} {data['RAs']}"
    dec = f"{data['DE-']}{data['DEd']} {data['DEm']} {data['DEs']}"
    return coord.SkyCoord(f'{ra} {dec}', unit=(u.hourangle, u.deg))

# Loop over all sources in the table
for source_data in source_table:
    # Make a SkyCoord object
    c = skycoord_from_table_row(source_data)
    # Perform a search around the specified coordinates
    r = requests.get(SST_URL,
                     params={**mipsgal_params, 'locstr': c.to_string()})

    # Extract the URL of the table of images
    img_tbl_url = xmltodict.parse(r.content)['result']['images']['metadata']
    # Need to switch to https and grab the file
    r2 = requests.get(img_tbl_url.replace('http:', 'https:'))
    # We need to remove the first line from the table so that it can be parsed
    table_lines = r2.content.decode().split('\n')[1:]
    # Then read it in as another astropy table
    img_table = Table.read(table_lines, format='ascii')

    # Select out all the images that are mosaic science images
    mosaic_images = {}
    for img_row in img_table:
        fname = img_row['fname']
        m_id = fname.split('/')[-1].split('_')[0]
        if img_row['file_type'] == 'science' and 'mosaics24' in fname:
            mosaic_images[m_id] = os.path.join(IMG_URL_ROOT, fname)

    # Now make postage stamps of all the selected images
    for m_id in mosaic_images:
        hdu = fits.open(mosaic_images[m_id])[0]
        w = WCS(hdu)
        # pixel coord of source
        i0, j0 = np.round(c.to_pixel(w))
        # find pixel limits of cut-out image window around source
        xpix_scale, ypix_scale = np.abs(w.wcs.cdelt)
        di = np.round(0.5*IMAGE_SIZE_DEGREES/xpix_scale)
        dj = np.round(0.5*IMAGE_SIZE_DEGREES/ypix_scale)
        win = slice(int(j0 - dj), int(j0 + dj)), slice(int(i0 - di), int(i0 + di))
        # Construct a new HDU
        hdr_win = w.slice(win).to_header()
        for k in 'Seq', 'Name', 'R0', 'PA':
            hdr_win[k] = source_data[k]
        hdr_win['MIPSGAL'] = m_id
        hdr_win['ORIGURL'] = mosaic_images[m_id]
        imwin = fits.PrimaryHDU(
            data=hdu.data[win],
            header=hdr_win)
        # Construct a suitable file name
        imid = f"{source_data['Seq']:04d}-{source_data['Name']}-{m_id}"
        imfn = os.path.join(OUTPUT_IMAGE_DIR, f'{imid}.fits')
        imwin.writeto(imfn, overwrite=True)
