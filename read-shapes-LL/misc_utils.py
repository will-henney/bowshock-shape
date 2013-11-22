import datetime
import getpass
import sys
import os
import json

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
        elif platform == "linux2":
            system = "linux"
        else:
            system = "unknown"
    return "{}@{}".format(user, system)


def fits_dirs():
    """
    Find the root directory for various collections of FITS files
    """
    dropbox_root = os.path.expanduser("~/Dropbox")
    user = who_am_i()
    if user.endswith("crya"):
        large_fits_dir = "/fs/nas11/other0/will/Orion"
    else:
        if user == "will@mac":
            large_fits_dir = "/Users/will/OrionTreasury"
        else:
            raise(EnvironmentError, "Unrecognised user@system - please add to misc_utils.py")

    if user.startswith("will"):
        dropbox_folder = "JorgeBowshocks"
    elif user.startswith("jorge") or user.startswith("angel"):
        dropbox_folder = "ProplydMIR"

    small_fits_dir = os.path.join(dropbox_root, dropbox_folder, "HST")
    return {
        "Bally": os.path.join(large_fits_dir, "BallyACS"), 
        "Robberto ACS": os.path.join(large_fits_dir, "acs"),
        "Robberto WFPC2": os.path.join(large_fits_dir, "wfc"),
        "PC": small_fits_dir,
        "PC mosaic": small_fits_dir,
        "WFC mosaic": small_fits_dir
    }

