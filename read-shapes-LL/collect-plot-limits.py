"""
Collect all the plot limits for a set of IMAGE-arcdata.json files

This is that we have a backup copy in case the files get overwritten

"""
from __future__ import print_function
import os
import json
import glob
import argparse
import misc_utils

parser = argparse.ArgumentParser(
    description="""Collect all the plot limits to plot-limits.json""")

parser.add_argument("--dirs", type=str,
                    default="j8oc*_wcs",
                    help="Pattern matching folders to be searched for SOURCE-arcdata.json files")

parser.add_argument("--debug", action="store_true",
                    help="Print out verbose debugging info")

cmd_args = parser.parse_args()

pattern = os.path.join(cmd_args.dirs, "*-arcdata.json")
collectfile = "image-plot-limits.json"
try:
    collect_db = json.load(open(collectfile))
except IOError:
    collect_db = {}

for dbfile in glob.glob(pattern):
    db = json.load(open(dbfile))
    for section in db.values():
        try: 
            name = section["extracted FITS file"]
            collect_db[name] = section["plot limits"]
        except KeyError:
            pass
            

misc_utils.update_json_file(collect_db, collectfile)

















