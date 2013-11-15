import json
import matplotlib.pyplot as plt

data = json.load(open("bally-offsets.json"))

# Choose a suitable layout of the sub-plots
for nrows, ncols in (1, 1), (1, 2), (2, 2), (2, 3), (3, 3), (3, 4), (3, 5):
    if len(data) < nrows*ncols:
        break 
else:
    raise ValueError, "Whoa!! Too many plots: " + str(len(data))

f = plt.figure()
# Additional graph of all the data
axx = f.add_subplot(nrows, ncols, nrows*ncols)
limits = [-2.0, 0.5, -1.25, 1.25]
for i, (field, d) in enumerate(data.items()):
    ax = f.add_subplot(nrows, ncols, i+1)
    # Plot offsets of individual sources
    ax.plot(d["dx"], d["dy"], "o", alpha=0.5)
    axx.plot(d["dx"], d["dy"], "o", alpha=0.5)
    # Plot trimean and iqr as a red cross
    ax.errorbar(d["avdx"][0], d["avdy"][0],
                xerr=d["avdx"][1], yerr=d["avdy"][1],
                c="r", lw=2.0, zorder=100, alpha=0.7)
    # Plot a +/- 0.5 arcsec cross at the origin
    ax.plot([-0.5, 0.5], [0.0, 0.0], "-k", lw=3.0, zorder=0)
    ax.plot([0.0, 0.0], [-0.5, 0.5], "-k", lw=3.0, zorder=0)
    # Label the outliers
    for dx, dy, s in zip(d["dx"], d["dy"], d["sources"]):
        if abs(dx - d["avdx"][0]) > d["avdx"][1] or abs(dy - d["avdy"][0]) > d["avdy"][1]:
            ax.text(dx, dy, s, fontsize=4, alpha=0.6)
    ax.axis("equal")
    ax.axis(limits)
    ax.grid()
    ax.set_xlabel(r"$\Delta\alpha$")
    ax.set_ylabel(r"$\Delta\delta$")
    ax.set_title(r"Bally Field {}, $N = {}$".format(field, len(d["sources"])))

axx.plot([-0.5, 0.5], [0.0, 0.0], "-k", lw=3.0, zorder=0)
axx.plot([0.0, 0.0], [-0.5, 0.5], "-k", lw=3.0, zorder=0)
axx.axis("equal")
axx.axis(limits)
axx.grid()
axx.set_xlabel(r"$\Delta\alpha$")
axx.set_ylabel(r"$\Delta\delta$")
axx.set_title(r"Bally All Fields")

f.set_size_inches(ncols*5, nrows*5)
f.tight_layout()
f.savefig("bally-offsets.pdf")
