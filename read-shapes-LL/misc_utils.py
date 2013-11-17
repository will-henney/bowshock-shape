import datetime
import getpass
import sys
import os
import json

def run_info():
    """
    Return info about this run:
    1. Name of program
    2. User who is running it
    3. Date and time it was run
    """
    return "{}: <{}> {}".format(
        os.path.basename(sys.argv[0]), 
        getpass.getuser(), 
        datetime.datetime.now().strftime('%Y-%m-%d %H:%M'),
    )


def update_json_file(data, jsonfile):
    """Try and update the JSON database file in the safest way possible

    This can still go wrong if the file is modified on disk in between
    the json.load and the json.dump

    """
    filedata = json.load(open(jsonfile))
    filedata.update(data)
    with open(jsonfile, "w") as f:
        json.dump(filedata, f, indent=4)


def short_image_name(long_name):
    if "hlsp_orion_hst_" in long_name:
        short_name = long_name.replace("hlsp_orion_hst_", 
                                       "Robberto_").replace("_v1_drz", "")
    elif long_name.startswith("j8oc"):
        short_name = long_name.replace("j8oc", "Bally_").replace("010_wcs", "")
    else:
        short_name = long_name
    return short_name
