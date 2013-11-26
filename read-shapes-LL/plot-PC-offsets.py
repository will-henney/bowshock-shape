import os
import glob
import matplotlib.pyplot as plt
import json

data = json.load(open("PC_offsets.json"))


#find a suitable layout for subplots

for nrows,ncols in (1,1),(1,2),(2,2),(2,3),(3,3),(3,4),(3,5):
    if len(data)< nrows*ncols:
        break
else:
    raise ValueError, "@#%#$, Too many plots =( " +str(len(data))

f =plt.figure()

axx= f.add_subplot(nrows,ncols,nrows*ncols)

for i,(field,d) in enumerate(data.items()):
    ax= f.add_subplot(nrows,ncols,i+1)
    #plot offsets for indivdual sources in the field subplot and the global subplot
    ax.plot(d["dx"],d["dy"],"o",alpha=0.5)
    axx.plot(d["dx"],d["dy"],"o",alpha=0.5)
    #Plot trimean and iqr as a red errorbar
    ax.errorbar(d["avdx"][0],d["avdy"][0],
                xerr=d["avdx"][1],yerr=d["avdy"][1],
                c="r",lw=2.0,zorder=100,alpha=0.7)
    # plot a black cross at origin
    ax.errorbar(0.0,0.0,xerr=0.5,yerr=0.5,c="k",lw=3.0,zorder=0.0)
    #Label the outliers
    for dx,dy,s in zip(d["dx"],d["dy"],d["sources"]):
        if abs(dx - d["avdx"][0]) > d["avdx"][1] or abs(dy - d["avdy"][0]) > d["avdy"][1]:
            ax.text(dx, dy, s, fontsize=4, alpha=0.6)
    ax.axis("equal")
    ax.grid()
    ax.set_xlabel(r"$\Delta\alpha$")
    ax.set_ylabel(r"$\Delta\delta$")
    ax.set_title(r"PC Field {}, $N = {}$".format(field, len(d["sources"])))

axx.errorbar(0.0,0.0,xerr=0.5,yerr=0.5,c="k",lw=3.0,zorder=0.0)
axx.axis("equal")
axx.grid()
axx.set_xlabel(r"$\Delta\alpha$")
axx.set_ylabel(r"$\Delta\delta$")
axx.set_title("PC all fields")

f.set_size_inches(ncols*5,nrows*5)
f.tight_layout()
f.savefig("PC-offsets.pdf")

