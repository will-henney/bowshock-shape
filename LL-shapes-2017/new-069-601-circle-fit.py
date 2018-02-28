import sys 
sys.path.append("/Users/will/Dropbox/circle-fit")
import circle_fit

for arc in "inner", "outer", "ridge":
    for delta_theta in 45, 50, 55, 60, 65, 70, 75, 80:
        plotfile = f"069-601-{arc}-{delta_theta:02d}.pdf"
        print('#### '*10)
        print("Creating", plotfile)
        circle_fit.plot_solution(
            f"new-069-601-{arc}.reg",
            "w069-601-Bally_01-extract.fits",
            plotfile,
            delta_theta=delta_theta, vmin=6.0, vmax=8.2, sigma=1.5)
