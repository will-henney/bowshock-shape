"""Find all the FITS images that contain each source

"""
from __future__ import print_function
import json
import os
import shutil
import glob
import argparse
import tempfile
import misc_utils
import fits_utils
import montage_wrapper
from astropy.io import fits
from astropy import wcs


def has_data(fitsfile, ra0, dec0):
    "Check whether there is any non-zero data in fitsfile at the source position (ra0, dec0)"
    fitspath = misc_utils.expand_fits_path(fitsfile)
    if cmd_args.debug: print("\nReading from", fitspath)
    hdu = fits_utils.get_image_hdu(fits.open(fitspath), debug=cmd_args.debug)
    if cmd_args.debug: print("Image size", hdu.data.shape)
    w = wcs.WCS(hdu.header)
    x, y = w.wcs_world2pix(ra0, dec0, 0)
    if cmd_args.debug: print("Coords", ra0, dec0, "->", x, y)
    i1 = int(x)
    j1 = int(y)
    pixel_data = hdu.data[j1:j1+2, i1:i1+2].mean()
    if cmd_args.debug: print("{}[{}:{},{}:{}] = {}".format(
            fitsfile, j1, j1+2, i1, i1+2, pixel_data))
    return pixel_data > 0.0 


def mcov_prefix():
    """Work around bug in some versions of mCoverageCheck"""
    if misc_utils.who_am_i() == "will@mac":
        return "BUG BUG "
    else:
        return ""


parser = argparse.ArgumentParser(
    description="""Find all the images that include a source""")
parser.add_argument("--dirs", type=str,
                    default="j8oc*_wcs",
                    help="Pattern matching folders to be searched for SOURCE-arcdata.json files")
parser.add_argument("--debug", action="store_true",
                    help="Print out verbose debugging info")
parser.add_argument("--full-check", action="store_true",
                    help="Check that each image has data of the source")

cmd_args = parser.parse_args()
pattern = os.path.join(cmd_args.dirs, "*-arcdata.json")

# Everything that is specific to a given user@machine is confined to
# the temporary directory.
tempdir = tempfile.mkdtemp()
#
# First, construct updated tables of all the FITS files available
#
tables = [os.path.join(tempdir, t) for t in ["small-all.tbl", "large-all.tbl"]]
for table, path in zip(tables, misc_utils.fits_dirs()):
    # This will rewrite the tables every time
    # FIXME should think of how to make it more efficient
    montage_wrapper.mImgtbl(path, table, recursive=True)

#
# Second, clean up and merge the tables
#
UNWANTED_DATASETS = ["wfi", "ispi", "spitzer"]
newtable_lines = []
for table in tables:
    table_lines = open(table).readlines()
    if not newtable_lines:
        # Copy header from first table
        newtable_lines = table_lines[0:3]
    for line in table_lines[3:]:
        fields = line.split()
        nhdu, fitsname = fields[-2:]
        # All the reasons to skip this line ...
        # ... we don't want ds9's backup files ...
        if ".bck.dir" in fitsname: continue
        # ... we only want the brightness, not the sigma or DQ ...
        if int(nhdu) > 1: continue
        # ... only some Bally images are corrected, skip the rest ...
        if "Bally" in fitsname and not "_wcs.fits" in fitsname: continue
        # ... some datasets we want to skip entirely ...
        for unwanted in UNWANTED_DATASETS:
            if "/" + unwanted + "/" in fitsname: continue
        # ... not sure what these are, but we don't want them ...
        if "/acs/" in fitsname and "_colorimage_" in fitsname: continue
        # ... try to select only the corrected images in the small dir ...
        if "small" in table:
            if not ("wcs_" in fitsname or "fix" in fitsname): continue
        # Otherwise, add to the combined list
        newtable_lines.append(line)

combo_table = os.path.join(tempdir, "combined-all.tbl")
with open(combo_table, "w") as f:
    f.writelines(newtable_lines)

#
# Third, loop through all the sources we can find, extracting the
# relevant images
#
imdb = {}
for dbfile in glob.glob(pattern):
    print(dbfile)
    # Load the source info database
    db = json.load(open(dbfile))
    # Extract relevant data
    source = db["star"]["id"]
    ra = db["star"]["RA_dg"]
    dec = db["star"]["Dec_dg"]

    # Look for source in both tables
    source_table = os.path.join(tempdir, "{}-images.tbl".format(source))
    if cmd_args.debug: print("\n\nWriting", source_table)
    # FIXME mCoverage check has a bug where it needs two extra args
    montage_wrapper.mCoverageCheck(mcov_prefix() + combo_table, source_table,
                                   "point", None, ra=ra, dec=dec)
    # Grab the names of the FITS images from the table that
    # mCoverageCheck wrote
    candidate_images = [misc_utils.contract_fits_path(line.split()[-1])
                        for line in open(source_table).readlines()
                        if ".fits" in line]
    # And optionally check that there is really data at the relevant position
    if cmd_args.full_check:
        imdb[source] = [fitsfile for fitsfile in candidate_images if has_data(fitsfile, ra, dec)]
    else:
        imdb[source] = candidate_images


#
# Fourth, write it out to a JSON file
#
misc_utils.update_json_file(imdb, "all-images.json")

# And clean up
if not cmd_args.debug:
    shutil.rmtree(tempdir)

# FIXME Check that each image has data at the position of the source

