import sys
import glob
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

try:
    ANGLE = sys.argv[1]
except:
    sys.exit(f"Usage: {sys.argv[0]} ANGLE")

zfiles = glob.glob(f"dust-bfield-stream-b{ANGLE}-z*.npz")

x, y = None, None
z, zplanes = [], []
for zfile in sorted(zfiles):
    with np.load(zfile) as data:
        zplanes.append(data['rho'])
        if x is None:
            x = data['x']
            y = data['y']
        z.append(int(zfile.split('.')[0][-4:])/100)
z = np.array(z)
# Add on the negative z planes by reflection symmetry about z=0
z = np.concatenate([-z[:0:-1], z])
zplanes = zplanes[:0:-1] + zplanes
cube = np.stack(zplanes)

xx = x[None, None, :]
yy = y[None, :, None]
zz = z[:, None, None]
rsq = xx**2 + yy**2 + zz**2


# Write out cube of density
fits.PrimaryHDU(
    data=cube
).writeto(f"dust-bfield-b{ANGLE}-cube.fits", overwrite=True)
# And cube of density * flux (Assuming ~ 1/R^2)
fits.PrimaryHDU(
    data=cube/rsq
).writeto(f"dust-bfield-b{ANGLE}-cube-F.fits", overwrite=True)
# And cube of density * T (Assuming T^4 ~ Flux)
fits.PrimaryHDU(
    data=cube/rsq**0.25
).writeto(f"dust-bfield-b{ANGLE}-cube-T.fits", overwrite=True)
# And cube of density * T^2 (Assuming T^4 ~ Flux)
fits.PrimaryHDU(
    data=cube/rsq**0.5
).writeto(f"dust-bfield-b{ANGLE}-cube-T2.fits", overwrite=True)
