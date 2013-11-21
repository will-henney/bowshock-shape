import json
import matplotlib.pyplot as plt

data = json.load(open("mosaic_offsets.json"))

#plot source offsets
plt.plot(data["dx"],data["dy"],"bo",alpha=0.5)

#plot trimean and iqr as a red cross

plt.errorbar(data["avdx"][0],data["avdy"][0],xerr=data["avdx"][1],
    yerr=data["avdy"][1],c="r",lw=2.0,zorder=100,alpha=0.7)

#Plot a cross at origin

plt.plot([-0.5, 0.5], [0.0, 0.0], "-k", lw=3.0, zorder=0)
plt.plot([0.0, 0.0], [-0.5, 0.5], "-k", lw=3.0, zorder=0)

for dx,dy,s in zip(data["dx"],data["dy"],data["sources"]):
    if abs(dx - data["avdx"][0]) > data["avdx"][1] or abs(dy - data["avdy"][0] > data["avdy"][1]):
        plt.text(dx,dy,s,fontsize=6,alpha=0.6)
plt.axis("equal")
plt.grid()
plt.xlabel(r"$\Delta\alpha$")
plt.ylabel(r"$\Delta\delta$")
plt.title("GO5469PCf* mosaic")
plt.savefig("mosaic-offsets.pdf")
