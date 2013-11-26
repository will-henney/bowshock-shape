from __future__ import print_function
import json
import argparse
import misc_utils

parser = argparse.ArgumentParser(
    description="""Remove absolute path references from JSON file""")

parser.add_argument("source", type=str,
                    default="LL1",
                    help="Name of source, taken as prefix for SOURCE-arcdata.json")
parser.add_argument("--user", type=str,
                    default="angel@crya",
                    help="User@system ID use for guessing paths")
parser.add_argument("--debug", action="store_true",
                    help="Print out verbose debugging info")

cmd_args = parser.parse_args()

dbfile = cmd_args.source + "-arcdata.json"

db = json.load(open(dbfile))

for secname, section in db.items():
    if "original FITS file" in section:
        oldpath = section["original FITS file"]
        newpath = misc_utils.contract_fits_path(oldpath, user=cmd_args.user)
        if cmd_args.debug:
            print("Contracting full paths in", secname)
            print("From", oldpath, "to", newpath)
            
        section["original FITS file"] = newpath

db["info"]["history"].append("Paths corrected by " + misc_utils.run_info())

misc_utils.update_json_file(db, dbfile)






