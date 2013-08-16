import numpy as np 
from equation6 import Shell
import matplotlib.pyplot as plt

N = 20                         # number of angles
theta = np.linspace(0.0, np.pi, N)

# Test isotropic and proplyd shells - separate figures
for innertype in "isotropic", "proplyd":
    print "**** {} ****".format(innertype)
    print
    for beta in 1.0, 0.2, 0.05, 0.02, 0.005:
        shell = Shell(beta=beta, innertype=innertype)
        R = shell.radius(theta)
        print R
        # set all non-positive values to NaN to prevent plotting
        R[R <= 0.0] = np.nan
        x, y = R*np.cos(theta), R*np.sin(theta)
        plt.plot(x, y, label="beta = {}".format(beta))

    plt.axis([-0.6, 0.6, -0.1, 1.1])
    plt.legend()
    plt.xlabel("z")
    plt.ylabel("r")
    plt.title("Shells from {} wind--isotropic wind interaction".format(innertype))
    plt.axes().set_aspect("equal")
    plt.savefig("shell-test-{}.png".format(innertype))
    plt.clf()
    print


# Now combine them on the same figure for small values of eta
N = 100                         # number of angles
theta = np.linspace(0.0, np.pi, N)
print
print "Combined figure"
for beta, color in zip([0.002, 0.005, 0.02], "rgb"):
    print "beta = ", beta
    shell = Shell(beta=beta, innertype="isotropic")
    R = shell.radius(theta)
    R[R <= 0.0] = np.nan
    x, y = R*np.cos(theta), R*np.sin(theta)
    plt.plot(x, y, color + "--", lw=2, label="isotropic: beta = {}".format(beta))

    shell = Shell(beta=beta, innertype="proplyd")
    R = shell.radius(theta)
    R[R <= 0.0] = np.nan
    x, y = R*np.cos(theta), R*np.sin(theta)
    plt.plot(x, y, color, lw=2, label="proplyd: beta = {}".format(beta))

plt.plot([-0.2, 0.25], [0.0, 0.0], "k")
plt.axis([-0.2, 0.25, -0.02, 0.35])
plt.legend(fontsize="small")
plt.xlabel("z")
plt.ylabel("r")
plt.title("Comparison between proplyd and isotropic inner winds".format(innertype))
plt.axes().set_aspect("equal")
plt.savefig("shell-test-compare.pdf")
plt.clf()


# Now make a graph of R0 vs R90
N = 10
# This time, we can take big steps and only go up to 90 degrees
theta = np.linspace(0.0, np.pi/2, N)
R0 = list()                     # the on-axis radius
R90 = list()                    # the perpendicular radius for isotropic inner wind
R90p = list()                    # the perpendicular radius for proplyd inner wind
beta = np.logspace(-3.0, -1.0)   # from 0.001 to 0.1 
print
print "R90 vs R0 figure"
for b in beta:
    # first do the isotropic wind
    shell = Shell(beta=b, innertype="isotropic")
    R = shell.radius(theta)
    R0.append(R[0])             # first element is R0
    R90.append(R[-1])           # last element is R90

    # second do the proplyd wind
    shell = Shell(beta=b, innertype="proplyd")
    R = shell.radius(theta)
    R90p.append(R[-1])           # last element is R90

 # Plotting the observationa data of LV1-LV5
R0LVs = np.array([0.16,0.24,0.29,0.17,0.21])
R90LVs = np.array([0.31,0.3,0.46,0.29,0.21])
LVs = ['LV1','LV2','LV3','LV4','LV5']
plt.axis([0,0.4,0,0.7])
plt.plot(R0LVs, R90LVs,'or')

#Plotting the theorical models
plt.plot(R0, R90, label="isotropic inner wind")
plt.plot(R0, R90p, label="proplyd inner wind")
#Labels of curves,axis and points
plt.legend()
plt.xlabel("R_0 / D")
plt.ylabel("R_90 / D")
for s in range(0,5):
    plt.figtext(2.5*R0LVs[s],10./7.*R90LVs[s],LVs[s],fontsize=10,ha='center')
plt.title("Perpendicular versus parallel bowshock radii")
plt.savefig("shell-test-R0-R90.pdf")

