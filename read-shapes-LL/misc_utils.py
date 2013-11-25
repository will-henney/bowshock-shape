import datetime
import getpass
import sys
import os
import json
import parse
import re

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
    imset, imset_data = find_image_set_re(long_name)
    return "_".join([imset.replace(" ", "_")] + [str(v) for v in imset_data.values ()])



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


IMAGE_SET_PATTERNS = {
    "Bally": "j8oc{field}010_wcs", 
    "Robberto ACS": "hlsp_orion_hst_acs_strip{field}_{filter}_v1_drz",
    "Robberto WFPC2": "hlsp_orion_hst_wfpc2_{field}_{filter}_v1_sci",
    "PC": "{field:d}f{filter:3w}",
    "PC mosaic": "wcs_GO5469PC{filter}e",
    "WFC mosaic": "mosaic{filter}",
}
UNKNOWN_SET_PATTERNS = {
    "Unknown dataset": "f{filter:3w}",
}


# Filters are optional "f" followed by 3 digits, followed by optional letters
RE_FILTER_NAME = "(?P<filter>f?\d\d\d[nwml]?p?)"
# Fields are 1 or 2 hexadecimal digits, or a digit followed by "l" or "r"
RE_FIELD_NAME = "(?P<field>[0-9a-flr]{1,2})"
IMAGE_SET_RE_PATTERNS = {
    "Bally": "j8oc" + RE_FIELD_NAME + "010_wcs", 
    "Robberto ACS": "hlsp_orion_hst_acs_strip" + RE_FIELD_NAME + "_" + RE_FILTER_NAME + "_v1_drz",
    "Robberto WFPC2": "hlsp_orion_hst_wfpc2_" +  RE_FIELD_NAME + "_" + RE_FILTER_NAME + "_v1_sci",
    "PC": RE_FIELD_NAME + RE_FILTER_NAME,
    "PC mosaic": "wcs_GO5469PC" + RE_FILTER_NAME + "e",
    "WFC mosaic": "mosaic" + RE_FILTER_NAME,
}
UNKNOWN_SET_RE_PATTERNS = {
    "Unknown dataset": RE_FILTER_NAME,
}




def find_image_set_re(filename):
    """Auto detect which dataset a given file comes from
    
    Returns name of image set and dict containing other info gleaned
    from the filename, such as field id and filter name.

    Rewrite using the re module
    """
    for imset, pattern in IMAGE_SET_RE_PATTERNS.items():
        match = re.match(pattern, filename)
        if match:
            return imset, match.groupdict()
    else:
        for imset, pattern in UNKNOWN_SET_RE_PATTERNS.items():
            match = re.search(pattern, filename)
            if match:
                return imset, match.groupdict()
        else:
            return "Unknown dataset no filter", {}


def find_image_set(filename):
    """Auto detect which dataset a given file comes from
    
    Returns name of image set and dict containing other info gleaned
    from the filename, such as field id and filter name.

    Uses the parse module, which is a little nicer than regexps,
    although maybe I should just learn to use named groups: 
    (?P<name>...), then accesed via match.groupdict()

    """
    for imset, pattern in IMAGE_SET_PATTERNS.items():
        match = parse.search(pattern, filename)
        if match:
            return imset, match.named
    else:
        for imset, pattern in UNKNOWN_SET_PATTERNS.items():
            match = parse.search(pattern, filename)
            if match:
                return imset, match.named
        else:
            return "Unknown dataset no filter", {}
    

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
            raise(EnvironmentError, 
                  "Unrecognised user@system - please add to misc_utils.py")

    if user.startswith("will"):
        dropbox_folder = "JorgeBowshocks"
    elif user.startswith("jorge") or user.startswith("angel"):
        dropbox_folder = "ProplydMIR"

    small_fits_dir = os.path.join(dropbox_root, dropbox_folder, "HST")
    return small_fits_dir, large_fits_dir

    # return {
    #     "Bally": os.path.join(large_fits_dir, "BallyACS"), 
    #     "Robberto ACS": os.path.join(large_fits_dir, "acs"),
    #     "Robberto WFPC2": os.path.join(large_fits_dir, "wfc"),
    #     "PC": small_fits_dir,
    #     "PC mosaic": small_fits_dir,
    #     "WFC mosaic": small_fits_dir
    # }

