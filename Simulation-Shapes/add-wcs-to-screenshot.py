import sys
from astropy.io import fits
from astropy.wcs import WCS

try: 
    prefix = sys.argv[1]
except:
    sys.exit(f"Usage: {sys.argv[0]} PREFIX")

channels = ['red', 'green', 'blue']
for channel in channels:
    hdulist = fits.open(f"{prefix}-{channel}.fits", mode="update")
    w = WCS(naxis=2)
    w.wcs.crpix = [0, 0]
    w.wcs.cdelt = [1.0/3600, 1.0/3600]
    w.wcs.crval = [30, -60]
    w.wcs.ctype = ["RA---AIR", "DEC--AIR"]
    hdulist[0].header = w.to_header()
    hdulist.flush()
    hdulist.close()
