import json

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

shelldata = json.load(open("rc-r0.json"))

params = {
    "font.family": "serif",
    "text.usetex": True,
    "text.latex.preamble": [r"\usepackage[varg]{txfonts}"],
    "figure.figsize": (5, 5),
}
matplotlib.rcParams.update(params)

lw = dict(isotropic=1, proplyd=2)
opacity = dict(isotropic=0.3, proplyd=0.7)

colors = "bgrmkcy"


for inn, modeldata in shelldata.items():
    for b, col in zip(modeldata, colors[:len(modeldata)]):
        inclinations = np.array(modeldata[b]["inc"])
        R0 = np.array(modeldata[b]["R0'"])
        Rc = np.array(modeldata[b]["Rc"])

        # Mask to select inclinations close to multiples of 15 degrees
        every15 = np.zeros(R0.shape, dtype=bool)
        for thisinc in 0.0, 15.0, 30.0, 45.0, 60.0, 75.0:
            iclosest = np.argmin(np.abs(np.degrees(inclinations) - thisinc))
            every15[iclosest] = True

        label = r'\(\beta={}\)'.format(b) if inn == "proplyd" else ""
        Y = Rc/R0
        # First, plot a line with all the inclinations
        plt.plot(R0, Y, '-', linewidth=lw[inn], c=col,
                 label=label, alpha=opacity[inn])
        # Second, add symbols every 15 degrees
        plt.plot(R0[every15], Y[every15], '.', c=col, alpha=opacity[inn])


# Add the observations to the plot
for pset, palpha in [
        ["best", 1.0],
        # ["second", 0.5]
]:
    pdata = np.genfromtxt(pset + "-proplyds.dat", names=True, dtype=None)
    plt.errorbar(pdata["R0D"], pdata["RcR0"],
                 xerr=pdata["ER0D"], yerr=pdata["ERcR0"],
                 fmt="ko", label="", alpha=palpha)
    for label, x, y, dx, dy in zip(
            pdata["Proplyd"], pdata["R0D"], pdata["RcR0"],
            pdata["dx"], pdata["dy"]):
        ha = "right" if dx < 0 else "left"
        va = "top" if dy < 0 else "bottom"
        plt.annotate(
            "{}".format(label),
            xy=(x, y), xytext=(dx, dy), fontsize="xx-small", alpha=palpha,
            textcoords = 'offset points', ha = ha, va = va,
            bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5*palpha),
            arrowprops = dict(arrowstyle='->', connectionstyle='arc3,rad=0')
        )


plt.xlabel(r"\(R'_0 / D'\)")
# avoid numbering the origin
epsilon = 1.e-6
plt.xlim(0.0 + epsilon, 0.5)
plt.ylim(0.0 - epsilon, 4 + epsilon)
plt.ylabel(r"\(R'_{c} / R'_0\)")
plt.legend(loc="best", ncol=2, prop=dict(size="x-small"))

plt.title("Curvature radii versus parallel bowshock radii")
plt.savefig("proplyd-shell-R0-Rc-errors.pdf")
plt.clf()
