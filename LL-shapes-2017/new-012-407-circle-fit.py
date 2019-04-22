import sys 
sys.path.append("/Users/will/Dropbox/circle-fit")
import circle_fit

id_ = "w012-407"
for arc in ["ridge"]:
    for delta_theta in 45, 50, 55, 60, 65, 70, 75, 80:
        plotfile = f"{id_}-{arc}-{delta_theta:02d}.pdf"
        print('#### '*10)
        print("Creating", plotfile)
        circle_fit.plot_solution(
            f"new-{id_}-{arc}.reg",
            f"{id_}-Bally_01-extract.fits",
            plotfile,
            delta_theta=delta_theta, vmin=3.0, vmax=4.3, sigma=5.0)
