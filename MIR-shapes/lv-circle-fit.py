import sys 
sys.path.append("/Users/will/Dropbox/circle-fit")
import circle_fit

fits_file = "/Users/will/Dropbox/JorgeBowshocks/MIR/SmithBally/wcs.mos11jy.fits"

try:
    source = str(sys.argv[1])
except:
    sys.exit(f"Usage: {sys.argv[0]} SOURCE")

for delta_theta in 65, 70, 75, 80:
    plotfile = f"{source}-{delta_theta:02d}.pdf"
    print('#### '*10)
    print("Creating", plotfile)
    circle_fit.plot_solution(
        f"{source}-smith-2005-forma.reg",
        fits_file,
        plotfile,
        delta_theta=delta_theta,
        vmin=0.0, vmax=0.01,
    )
