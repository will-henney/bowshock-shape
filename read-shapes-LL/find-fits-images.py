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
newtable_lines = []
for table in tables:
    table_lines = open(table).readlines()
    if not newtable_lines:
        # Copy header from first table
        newtable_lines = table_lines[0:3]
    for line in table_lines[3:]:
        fields = line.split()
        nhdu, fitsname = fields[-2:]
        # All the reasons to skip this line
        if ".bck.dir" in fitsname: continue
        if int(nhdu) > 1: continue
        if "Bally" in fitsname and not "_wcs.fits" in fitsname: continue
        if "/wfi/" in fitsname or "/ispi/" in fitsname: continue
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
    # Load the source info database
    db = json.load(open(dbfile))
    # Extract relevant data
    source = db["star"]["id"]
    ra = db["star"]["RA_dg"]
    dec = db["star"]["Dec_dg"]

    # Look for source in both tables
    source_table = os.path.join(tempdir, "{}-images.tbl".format(source))
    if cmd_args.debug: print("Writing", source_table)
    # FIXME mCoverage check has a bug where it needs two extra args
    montage_wrapper.mCoverageCheck("BUG BUG " + combo_table, source_table,
                                   "point", None, ra=ra, dec=dec)
    # Grab the names of the FITS images from the table that
    # mCoverageCheck wrote
    imdb[source] = [misc_utils.contract_fits_path(line.split()[-1])
                    for line in open(source_table).readlines()
                    if ".fits" in line]



#
# Fourth, write it out to a JSON file
#
misc_utils.update_json_file(imdb, "all-images.json")

# And clean up
if not cmd_args.debug:
    shutil.rmtree(tempdir)

# FIXME Check that each image has data at the position of the source

