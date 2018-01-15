import sys
import glob
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

try:
    PREFIX = sys.argv[1]
except:
    sys.exit(f"Usage: {sys.argv[0]} PREFIX")

zfile = PREFIX + ".npz"

with np.load(zfile) as data:
    cube = data['rho'].T
    x = data['x']
    y = data['y']
    z = data['z']

xx = x[None, None, :]
yy = y[None, :, None]
zz = z[:, None, None]
rsq = xx**2 + yy**2 + zz**2

w = WCS(naxis=3)
w.wcs.crpix = [1, 1, 1]
w.wcs.crval = [x[0], y[0], z[0]]
w.wcs.cdelt = [x[1] - x[0], y[1] - y[0], z[1] - z[0]]

# Write out cube of density
fits.PrimaryHDU(
    data=cube, header=w.to_header(),
).writeto(f"{PREFIX}-cube.fits", overwrite=True)
# And cube of density * flux (Assuming ~ 1/R^2)
fits.PrimaryHDU(
    data=cube/rsq, header=w.to_header(),
).writeto(f"{PREFIX}-cube-F.fits", overwrite=True)
# And cube of density * T (Assuming T^4 ~ Flux)
fits.PrimaryHDU(
    data=cube/rsq**0.25, header=w.to_header(),
).writeto(f"{PREFIX}-cube-T.fits", overwrite=True)
# And cube of density * T^2 (Assuming T^4 ~ Flux)
fits.PrimaryHDU(
    data=cube/rsq**0.5, header=w.to_header(),
).writeto(f"{PREFIX}-cube-T2.fits", overwrite=True)
