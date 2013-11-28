"""Find all the FITS images that contain each source

"""
from __future__ import print_function
import json
import os
import glob
import argparse
import misc_utils
import montage_wrapper


parser = argparse.ArgumentParser(
    description="""Find all the images that include a source""")
parser.add_argument("--dirs", type=str,
                    default="j8oc*_wcs",
                    help="Pattern matching folders to be searched for SOURCE-arcdata.json files")
parser.add_argument("--debug", action="store_true",
                    help="Print out verbose debugging info")

cmd_args = parser.parse_args()
pattern = os.path.join(cmd_args.dirs, "*-arcdata.json")

#
# First, construct updated tables of all the FITS files available
#
tables = "small-all.tbl", "large-all.tbl"
for table, path in zip(tables, misc_utils.fits_dirs()):
    # This will rewrite the tables every time
    # FIXME should think of how to make it more efficient
    montage_wrapper.mImgtbl(path, table, recursive=True)

#
# Second, loop through all the sources we can find, extracting the
# relevant images
#
for dbfile in glob.glob(pattern):
    # Load the source info database
    db = json.load(open(dbfile))
    # Extract relevant data
    source = db["star"]["id"]
    ra = db["star"]["RA_dg"]
    dec = db["star"]["Dec_dg"]

    # Look for source in both tables
    for table in tables:
        source_table = table.replace("-all.tbl", "-{}.tbl".format(source))
        # FIXME mCoverage check has a bug where it needs two extra args
        montage_wrapper.mCoverageCheck("BUG BUG " + table, source_table,
                                       "point", None, ra=ra, dec=dec)


# Filter out the images we want

# Check that each image has data at the position of the source

