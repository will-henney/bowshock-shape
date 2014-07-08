from __future__ import print_function
import datetime
import getpass
import sys
import os
import json
import re
import numpy as np

def run_info():
    """
    Return info about this run:
    1. Name of program
    2. User@system who is running it
    3. Date and time it was run
    """
    return "{}: <{}> {}".format(
        os.path.basename(sys.argv[0]), 
        who_am_i(), 
        datetime.datetime.now().strftime('%Y-%m-%d %H:%M'),
    )


class NumpyAwareJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray) and obj.ndim <= 1:
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def update_json_file(data, jsonfile):
    """Try and update the JSON database file in the safest way possible

    This can still go wrong if the file is modified on disk in between
    the json.load and the json.dump

    """
    try:
        filedata = json.load(open(jsonfile))
    except IOError:
        filedata = {}
    filedata.update(data)
    with open(jsonfile, "w") as f:
        json.dump(filedata, f, indent=4, cls=NumpyAwareJSONEncoder)


def short_image_name(long_name):
    imset, imset_data = find_image_set(long_name)
    return "_".join([imset.replace(" ", "_")] + [str(v) for v in imset_data])


def who_am_i():
    """
    Find who I am: user@system

    This is used to find where some of the data files are
    """
    user = getpass.getuser()
    if os.path.exists("/fs/nas11/other0/will"):
        # If we can see the NAS disk we don't care which OS
        system = "crya"
    else:
        # If not, we are probably on someone's laptop
        platform = os.sys.platform
        if platform == "darwin":
            system = "mac"
        elif "linux" in platform:
            system = "linux"
        else:
            system = "unknown"
    return "{}@{}".format(user, system)


# Filters are "f" followed by 3 digits, followed by optional letters
FILTER = "(?P<filter>f\d\d\d[nwml]?p?)"
# Fields are 1 or 2 hexadecimal digits, or a digit followed by "l" or "r"
FIELD = "(?P<field>[0-9a-flr]{1,2})"
IMAGE_SET_PATTERNS = {
    "Bally": "j8oc" + FIELD + "010_wcs", 
    "Robberto ACS": "hlsp_orion_hst_acs_strip" + FIELD + "_" + FILTER + "_v1_drz",
    "Robberto WFPC2": "hlsp_orion_hst_wfpc2_" +  FIELD + "_" + FILTER + "_v1_sci",
    "PC": "wcs_" + FIELD + FILTER,
    "PC mosaic": "wcs_GO5469PC" + FILTER + "e",
    "WFC mosaic": "mosaic" + FILTER,
    "ACS Ramp fr505n": "oiii-trap-fix",
}
UNKNOWN_SET_PATTERNS = {
    "Unknown dataset": FILTER,
}


def find_image_set(filename):
    """Auto detect which dataset a given file comes from
    
    Returns name of image set and list containing other info gleaned
    from the filename, such as field id and filter name.

    Rewrite using the re module - wasn't too hard
    """
    for imset, pattern in IMAGE_SET_PATTERNS.items():
        match = re.match(pattern, filename)
        if match:
            return imset, match.groups()
    else:
        for imset, pattern in UNKNOWN_SET_PATTERNS.items():
            match = re.search(pattern, filename)
            if match:
                return imset, match.groups()
        else:
            return "Unknown dataset no filter", []

#   
# Utility functions to deal with data directories
#
def fits_dirs(user=None):
    """
    Find the root directory for various collections of FITS files
    """
    dropbox_root = os.path.expanduser("~/Dropbox")
    if user is None:
        user = who_am_i()
    if user.endswith("crya"):
        large_fits_dir = "/fs/nas11/other0/will/Orion"
    else:
        if user == "will@mac":
            large_fits_dir = "/Users/will/Work/OrionTreasury"
        elif user.endswith("linux"):
            print("No large fits directory") #Avoid Error when using Laptop
            large_fits_dir= ""
        else:
            raise EnvironmentError("Unrecognised user@system - please add to misc_utils.py")
    if user.startswith("will"):
        dropbox_folder = "JorgeBowshocks"
    elif user.split("@")[0] in ["jtarango", "dark_lxndrs"]:
        dropbox_folder = "ProplydMIR"
    elif user.split("@")[0] in ["angel", "luis"]:
        dropbox_folder = "LuisBowshocks"

    
    # Special cases for Dropbox location on CRyA machines
    if user == "jtarango@crya":
        dropbox_root ="/fs/posgrado01/other0/jtarango/Dropbox"
    if user == "angel@crya":
        dropbox_root = "/fs/posgrado01/other0/angel/Dropbox" 

    small_fits_dir = os.path.join(dropbox_root, dropbox_folder, "HST")
    return small_fits_dir, large_fits_dir


def expand_fits_path(p, user=None):
    """
    Expand out to the full path
    """
    small_fits_dir, large_fits_dir = fits_dirs(user)
    if "$SMALL" in p:
        return p.replace("$SMALL_FITS_DIR", small_fits_dir)
    elif "$LARGE" in p:
        return p.replace("$LARGE_FITS_DIR", large_fits_dir)
    else:
        return p


def contract_fits_path(p, user=None):
    """
    Canonicalize path
    """
    small_fits_dir, large_fits_dir = fits_dirs(user)
    return p.replace(small_fits_dir, "$SMALL_FITS_DIR").replace(large_fits_dir, "$LARGE_FITS_DIR")


def image_set_fits_dir(image_name):
    """Find the full path a given set of image files

    We are quite flexible with the format of image_name.  It could be
    a short image name 

    """
    small_fits_dir, large_fits_dir = fits_dirs()
    if "Bally" in image_name:
        return os.path.join(large_fits_dir, "Bally-ACS")
    elif "Robberto" in image_name:
        if "ACS" in image_name:
            return os.path.join(large_fits_dir, "acs")
        elif "WFPC2" in image_name:
            return os.path.join(large_fits_dir, "wfpc2")
        else:
            raise NotImplementedError()
    else:
        return small_fits_dir
