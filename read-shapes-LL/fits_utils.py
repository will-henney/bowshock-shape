from astropy.io import fits

def get_image_hdu(hdulist, debug=False):
    """
    Return the main image from a FITS file

    Copy in the keywords from the Primary HDU if "INHERIT" = T
    """

    if len(hdulist) == 1:
        # If there is only one HDU in the FITS file then use it
        hdu, = hdulist
        if debug:
            print "Only one HDU found"
    else:
        # Otherwise, use the HDU named "SCI"
        hdu = hdulist["SCI"]
        if debug:
            print "Using 'SCI' HDU"
        # And copy in the keywords from the Primary HDU 
        if hdu.header.get("INHERIT"):
            hdu.header.update(hdulist[0].header)
            if debug:
                print "Inheriting from Primary HDU"
    return hdu



def get_instrument_configuration(hdu, debug=False):
    """
    Return a string that describes the instrument and camera
    """
    filter_kwds = ["FILTER1", "FILTNAM1"]
    camera = hdu.header.get("INSTRUME", "Unknown Camera")
    for k in filter_kwds:
        filtro = hdu.header.get(k)
        if filtro:
            break
    else:
        filtro = "Unknown Filter"
    return {"filter": filtro, "camera": camera}
    
