import sys 
sys.path.append("/Users/will/Dropbox/circle-fit")
import circle_fit

id_ = "LL4"
for arc in ["ridge"]:
    for delta_theta in 45, 50, 55, 60, 65, 70, 75, 80:
        plotfile = f"{id_}-{arc}-{delta_theta:02d}.pdf"
        print('#### '*10)
        print("Creating", plotfile)
        circle_fit.plot_solution(
            f"new-{id_}-{arc}.reg",
            f"{id_}-Bally_24-extract.fits",
            plotfile,
            delta_theta=delta_theta, vmin=0.85, vmax=1.0, sigma=5.0, maxiter=5)
