import numpy as np
import json
from astropy.io import fits
from astropy import wcs
from astropy import coordinates as coord
from astropy import units as u 
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from image_statistics import trimean_and_iqr, robust_statistics
from misc_utils import run_info, update_json_file
import fits_utils
import pyregion

def offset_circle_radius(theta, b):
    """
    Radius from origin to a unit circle that is offset a distance b (< 1) away from origin

    theta is angle from axis (direction opposite to b)
    """
    assert b < 1.0 and b >= 0.0, "Parameter b (= {:.3f}) must be in range [0, 1]".format(b)
    return np.sqrt(1 - b**2*np.sin(theta)**2) - b*np.cos(theta)




parser = argparse.ArgumentParser(
    description="""Calculate surface brightness profile for bowshock arcs""")

parser.add_argument("source", type=str,
                    default="LL1",
                    help="Name of source (prefix for files)")
parser.add_argument("--debug", action="store_true",
                    help="Print out verbose debugging info")

cmd_args = parser.parse_args()

dbfile = cmd_args.source + "-arcdata.json"
arcdata = json.load(open(dbfile))


for image_name in arcdata:
    if not "extracted FITS file" in arcdata[image_name]:
        # Jump over sections that are not images
        continue
    print "**************************************************"
    print "Calculating shell brightness for ", image_name
    print "**************************************************"

    fitsfile = arcdata[image_name]["extracted FITS file"]
    hdulist = fits.open(fitsfile)
    hdu = fits_utils.get_image_hdu(hdulist, debug=cmd_args.debug)

    w = wcs.WCS(hdu.header)
    # Pixel coordinate arrays
    x, y = np.meshgrid(np.arange(hdu.data.shape[1]), np.arange(hdu.data.shape[0]))
    # Sky coordinate arrays
    ra, dec = w.wcs_pix2world(x, y, 0)
    # Find coordinates of the central star
    ra0 = coord.Longitude(arcdata["star"]["RA"], unit=u.hour)
    dec0 = coord.Latitude(arcdata["star"]["Dec"], unit=u.deg)
    # Calculate sky coordinate offset arrays in arcsec with respect to star.
    # (Note that this will not work near the poles, or over a field that
    # spans a large angle of sky, but that does not bother us.)
    X = np.cos(np.radians(dec))*(ra*u.deg - ra0).to(u.arcsec)
    Y = (dec*u.deg - dec0).to(u.arcsec)
    # Radial separation, R, and position angle of separation, PA.
    R = np.hypot(X, Y)
    PA = coord.Longitude(np.arctan2(X, Y), unit=u.radian).to(u.degree)

    # Get circle parameters for outer/inner arc
    theta = {}
    R_arc = {}
    for arc in "inner", "outer":
        if arc in arcdata:
            PAc = coord.Longitude(arcdata[arc]["PAc"], u.deg)
            Rc = arcdata[arc]["Rc"] * u.arcsec
            Xc = arcdata[arc]["xc"] * u.arcsec
            Yc = arcdata[arc]["yc"] * u.arcsec
            b = np.hypot(Xc, Yc)/Rc
            # Angle from axis of arc
            # (in range [-180:180], measured anticlockwise)
            theta[arc] = (PA - PAc).wrap_at(180*u.deg)
            # For each pixel, this is the radius from star of the outer arc
            # (in circle approximation) for the pixel's theta
            R_arc[arc] = Rc*offset_circle_radius(theta[arc], b)
        else:
            raise NotImplementedError("Case of only one arc not covered yet.")

    # Position with respect to the arcs:
    # z < 0 -- inside inner arc
    # z = 0 -> 1 -- between arcs
    # z > 0 -- outside outer arc
    z = (R - R_arc["inner"]) / (R_arc["outer"] - R_arc["inner"])



    # Various masks for the shell, background, etc
    mask = {
        "core": R < 2.5*R_arc["outer"],
        "shell": np.abs(z - 0.5) <= 0.5,
        "shell center": np.abs(z - 0.5) <= 0.25,
        "<45": np.abs(theta["outer"].deg) <= 45,
        "<60": np.abs(theta["outer"].deg) <= 60,
        "<90": np.abs(theta["outer"].deg) <= 90,
        "<15":  np.abs(theta["outer"].deg) <= 15.0, 
        "+15 to +45":  np.abs(theta["outer"].deg - 30.0) <= 15.0, 
        "-15 to -45":  np.abs(theta["outer"].deg + 30.0) <= 15.0, 
        "bg": np.abs(z - 1.5) <= 0.3,
    }

    try:
        regions = pyregion.open(cmd_args.source + "-mask.reg")
    except IOError:
        print "No bad pixel mask found"
        mask["good"] = np.ones_like(x, dtype=bool)
    else:
        print "Using " + cmd_args.source + "-mask.reg as bad pixel mask"
        mask["good"] = ~regions.get_mask(hdu=hdu)
    #
    # Calculate robust statistics 
    avbg, wbg = trimean_and_iqr(hdu.data[mask["<45"] & mask["bg"] & mask["good"]])
    avsh, wsh = trimean_and_iqr(hdu.data[mask["<45"] & mask["shell"] & mask["good"]])
    avshc, wshc = trimean_and_iqr(hdu.data[mask["<45"] & mask["shell center"] & mask["good"]])
    ymin = avbg - 2*wbg 
    ymax = avsh + 2*wsh


    print "BG trimean = {:.2f}, iqr = {:.2f}".format(avbg, wbg)
    print "Shell trimean = {:.2f}, iqr = {:.2f}".format(avsh, wsh)
    print "Adopting plot range of {:.2f} to {:.2f}".format(ymin, ymax)

    # In the future, we will generalise this to do other filters/cameras 


    # Save brightness statistics to a new JSON file
    arcdata[image_name]["background"] = {"value": avbg, "delta": wbg}
    arcdata[image_name]["shell"] = {"value": avsh, "delta": wsh}
    arcdata[image_name]["shell center"] = {"value": avshc, "delta": wshc}
    arcdata["info"]["history"].append("Shell data for " + image_name + " added/modified by " + run_info())

    # Save brightness statistics in the help section
    help_brightness = {"shell":"Trimean brightness and interquartile (difference between quartiles) of shell",
                       "shell center":"Trimean brightness and interquartile (difference between quartiles) of shell center",
                       "background":"Trimean brightness and interquartile (difference between quartiles) of background"}

    arcdata["help"].update(brightness=help_brightness)

    update_json_file(arcdata, dbfile)

    # Calculate average binned profiles versus theta
    th_edges = np.linspace(-60.0, 60.0, 25)
    th_centers = 0.5*(th_edges[:-1] + th_edges[1:])
    m = mask["<60"] & mask["good"] & mask["bg"]
    bg, dbg = robust_statistics(theta["outer"][m].deg, hdu.data[m], th_edges)
    m = mask["<60"] & mask["good"] & mask["shell center"]
    sh, dsh = robust_statistics(theta["outer"][m].deg, hdu.data[m], th_edges)

    # Calculate average binned profiles versus radius
    z_edges = np.linspace(-2.0, 3.0, 50)
    z_centers = 0.5*(z_edges[:-1] + z_edges[1:])
    m = mask["good"] & mask["<15"]
    axis, daxis = robust_statistics(z[m], hdu.data[m], z_edges)
    m = mask["good"] & mask["+15 to +45"]
    upper, dupper = robust_statistics(z[m], hdu.data[m], z_edges)
    m = mask["good"] & mask["-15 to -45"]
    lower, dlower = robust_statistics(z[m], hdu.data[m], z_edges)

    # Plot graph of radial profiles
    #
    plot_prefix = "-".join([cmd_args.source, image_name, "arcbright"])
    plt.plot(z_centers, lower, color="y", alpha=0.8, zorder=111, lw=2.0, label=r"$\theta = -45^\circ$ to $-15^\circ$")
    plt.fill_between(z_centers, lower-dlower, lower+dlower, color="y", alpha=0.05, lw=0, zorder=101)
    plt.plot(z_centers, axis, color="c", alpha=0.8, zorder=112, lw=2.0, label=r"$\theta = -15^\circ$ to $+15^\circ$")
    plt.fill_between(z_centers, axis-daxis, axis+daxis, color="c", alpha=0.05, lw=0, zorder=102)
    plt.plot(z_centers, upper, color="b", alpha=0.8, zorder=113, lw=2.0, label=r"$\theta = +15^\circ$ to $+45^\circ$")
    plt.fill_between(z_centers, upper-dupper, upper+dupper, color="b", alpha=0.05, lw=0, zorder=103)

    m = mask["<60"] & mask["good"] & mask["core"]
    plt.scatter(z[m], hdu.data[m], s=2, c=theta["outer"][m].deg,
                marker="o", cmap=plt.cm.gist_rainbow, alpha=0.6, faceted=False, zorder=100)
    cb = plt.colorbar()
    cb.set_label("theta")
    plt.xlabel("Radius relative to shell: (R - R_in) / (R_out - R_in)")
    plt.ylabel("Surface brightness")
    plt.legend(prop={"size": "small"})
    plt.title(" ".join([cmd_args.source, arcdata[image_name]["camera"], arcdata[image_name]["filter"]]))
    plt.xlim(-2.0, 3.0)
    plt.ylim(ymin, ymax)
    plt.grid()
    plt.savefig(plot_prefix + "-z.png", dpi=600)


    #
    # Plot graph of angular profiles
    #
    plt.clf()
    plt.plot(th_centers, bg, color="k", alpha=0.3, zorder=110, label="Background")
    plt.fill_between(th_centers, bg-dbg, bg+dbg, color="k", alpha=0.1, lw=0, zorder=100)
    plt.plot(th_centers, sh, color="g", alpha=0.3, zorder=111, lw=2, label="Shell center")
    plt.fill_between(th_centers, sh-dsh, sh+dsh, color="m", alpha=0.1, lw=0, zorder=101)
    m = mask["<60"] & mask["good"] & mask["shell"]
    mb = mask["<60"] & mask["good"] & mask["bg"]
    plt.scatter(theta["outer"][m].deg, hdu.data[m], s=4, c=z[m], marker="o",
                cmap=plt.cm.gist_rainbow, alpha=0.6, faceted=False, zorder=20)
    cb = plt.colorbar()
    cb.set_label("z = (R - R_in) / (R_out - R_in)")
    plt.scatter(theta["outer"][mb].deg, hdu.data[mb], s=2, c="k", marker="o",
                cmap=plt.cm.gist_rainbow, alpha=0.2, faceted=False, zorder=10)
    plt.xlabel("Angle, theta, from outer shell axis")
    plt.ylabel("Surface brightness")
    plt.legend()
    plt.title(" ".join([cmd_args.source, arcdata[image_name]["camera"], arcdata[image_name]["filter"]]))
    plt.ylim(ymin, ymax)
    plt.grid()
    plt.savefig(plot_prefix + "-th.png", dpi=600)


