from __future__ import print_function

def get_image_hdu(hdulist, debug=False, fix=True):
    """
    Return the main image from a FITS file

    Copy in the keywords from the Primary HDU if "INHERIT" = T
    """

    if len(hdulist) == 1:
        # If there is only one HDU in the FITS file then use it
        hdu, = hdulist
        if debug:
            print("Only one HDU found")
    else:
        # Otherwise, use the HDU named "SCI"
        hdu = hdulist["SCI"]
        if debug:
            print("Using 'SCI' HDU")
        # And copy in the keywords from the Primary HDU 
        if hdu.header.get("INHERIT"):
            hdu.header.update(hdulist[0].header.items())
            if debug:
                print("Inheriting from Primary HDU")

    # Fix up potential problems with the header
    if fix:
        if "EQUINOX" in hdu.header:
            pass
    return hdu




def get_instrument_configuration(hdu, debug=False):
    """
    Return a string that describes the instrument and camera
    """
    filter_kwds = ["FILTER1", "FILTNAM1", "FILTER2", "FILTNAM2"]
    camera = hdu.header.get("INSTRUME", "Unknown Camera")
    date = hdu.header.get("DATE-OBS", "Unknown Date")
    for k in filter_kwds:
        filter_ = hdu.header.get(k)
        if filter_ and str(filter_).startswith("F"):
            break
    else:
        filter_ = "Unknown Filter"
    return {"filter": filter_, "camera": camera, "date": date}
    
