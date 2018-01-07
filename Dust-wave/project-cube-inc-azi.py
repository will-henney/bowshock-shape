import sys
import numpy as np
from scipy.interpolate import interpn
from astropy.io import fits
from astropy.wcs import WCS

try:
    PREFIX = sys.argv[1]
    INC = float(sys.argv[2])
    AZI = float(sys.argv[3])
except:
    sys.exit(f"Usage: {sys.argv[0]} PREFIX INC AZI")


# radians versions of inclination and azimuth
inc = np.radians(INC)
azi = np.radians(AZI)

#
# All vectors have components in body frame cartesian basis (i, j, k)
#

# Unit vectors for body frame (note that coordinates are [z, y, x],
# which is why we reverse them)
#
ivec = np.array([1.0, 0.0, 0.0])[::-1]
jvec = np.array([0.0, 1.0, 0.0])[::-1]
kvec = np.array([0.0, 0.0, 1.0])[::-1]

# Unit vectors for observer frame (primed) - again, all reversed to [z, y, x]
iprime = np.array([np.cos(inc),
                   0.0,
                   np.sin(inc)])[::-1]
jprime = np.array([-np.sin(inc)*np.sin(azi),
                   np.cos(azi),
                   np.cos(inc)*np.sin(azi)])[::-1]
kprime = np.array([-np.sin(inc)*np.cos(azi),
                   -np.sin(azi),
                   np.cos(inc)*np.cos(azi)])[::-1]


smin, smax = -5.0, 5.0

N = 201
# Coordinate grid in observer frame cartesian basis (i', j', k')
xprime = np.linspace(smin, smax, N)[None, None, :]
yprime = np.linspace(smin, smax, N)[None, :, None]
zprime = np.linspace(smin, smax, N)[:, None, None]

# Body frame coordinates of observer grid. Shape: (NZ, NY, NX, 3)
body_sample_coords = (
    xprime[:, :, :, None]*iprime[None, None, None, :]
    + yprime[:, :, :, None]*jprime[None, None, None, :]
    + zprime[:, :, :, None]*kprime[None, None, None, :])

# Read in body frame data
hdu, = fits.open(PREFIX + ".fits")
w = WCS(hdu.header)
nz, ny, nx = hdu.data.shape
xpts, _, _ = w.wcs_pix2world(range(nx), [0]*nx, [0]*nx, 0)
_, ypts, _ = w.wcs_pix2world([0]*ny, range(ny), [0]*ny, 0)
_, _, zpts = w.wcs_pix2world([0]*nz, [0]*nz, range(nz), 0)

data_samples = interpn(
    (zpts, ypts, xpts), hdu.data, body_sample_coords,
    bounds_error=False, fill_value=0.0)

image = data_samples.sum(axis=0)
wim = WCS(naxis=2)
wim.wcs.crpix = [1, 1]
wim.wcs.crval = [xprime[0, 0, 0], yprime[0, 0, 0]]
wim.wcs.cdelt = [xprime[0, 0, 1] - xprime[0, 0, 0],
                 yprime[0, 1, 0] - yprime[0, 0, 0]]

outfile = f"{PREFIX}-{int(INC):03d}-{int(AZI):03d}.fits"
fits.PrimaryHDU(
    data=image, header=wim.to_header()
).writeto(outfile, overwrite=True)

print(outfile, end="")
