import sys 
sys.path.append("/Users/will/Dropbox/circle-fit")
import circle_fit

for arc in "inner", "outer", "ridge":
    for delta_theta in 45, 50, 55, 60, 65, 70, 75, 80:
        plotfile = f"000-400-{arc}-{delta_theta:02d}.pdf"
        print('#### '*10)
        print("Creating", plotfile)
        circle_fit.plot_solution(
            f"new-w000-400-{arc}.reg",
            "w000-400-Bally_09-extract.fits",
            plotfile,
            delta_theta=delta_theta)
