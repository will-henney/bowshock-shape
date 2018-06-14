import sys 
sys.path.append("/Users/will/Dropbox/circle-fit")
import circle_fit

datapaths = [
    "/Users/will/Dropbox/OrionHST-2012/Combine/",
    "/Users/will/Dropbox/OrionBally-2016/data-historic/",
    "/Users/will/Dropbox/OrionBally-2016/data-2016-01-11/",
]

fits_files = {
    "oiii-1993-pc": datapaths[0] + "f502n-allpc-align-rob_drz_sci.fits",
    "oiii-2004-acs": datapaths[0] + "fr505n-5007-align-rob_drz_sci.fits", 
    "ha-1993-pc": datapaths[1] + "wcs_3f656.fits",
    "ha-2004-acs": datapaths[1] + "j8oc01010_wcs.fits",
    "ha-2005-acs": datapaths[1] + "hlsp_orion_hst_acs_strip0l_f658n_v1_drz.fits",
    "ha-2015-wfc3": datapaths[2] + "icaz01040_drz.fits"
    }

vlim = {
    "oiii-1993-pc": [0.24, 1.1],
    "oiii-2004-acs": [20, 60],
    "ha-1993-pc": [75, 118],
    "ha-2004-acs": [50, 75],
    "ha-2005-acs": [50, 75],
    "ha-2015-wfc3": [12.5, 19.3],
}

try:
    image_id = str(sys.argv[1])
except:
    sys.exit(f"Usage: {sys.argv[0]} SOURCE")

for delta_theta in 40, 45, 50, 55, 60, 65, 70, 75:
    plotfile = f"lv4-{image_id}-{delta_theta:03d}.pdf"
    print('#### '*10)
    print("Creating", plotfile)
    circle_fit.plot_solution(
        f"lv4-{image_id}-forma.reg",
        fits_files[image_id],
        plotfile,
        delta_theta=delta_theta,
        vmin=vlim[image_id][0], vmax=vlim[image_id][1],
        verbose=False, maxiter=5
    )
