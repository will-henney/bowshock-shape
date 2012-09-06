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


# Now combine them on the same figure
N = 100                         # number of angles
theta = np.linspace(0.0, np.pi, N)
for beta, color in zip([0.002, 0.005, 0.02], "rgb"):
    shell = Shell(beta=beta, innertype="isotropic")
    R = shell.radius(theta)
    R[R <= 0.0] = np.nan
    x, y = R*np.cos(theta), R*np.sin(theta)
    plt.plot(x, y, color + "--", label="isotropic: beta = {}".format(beta))

    shell = Shell(beta=beta, innertype="proplyd")
    R = shell.radius(theta)
    R[R <= 0.0] = np.nan
    x, y = R*np.cos(theta), R*np.sin(theta)
    plt.plot(x, y, color, label="proplyd: beta = {}".format(beta))

plt.axis([-0.2, 0.2, -0.05, 0.35])
plt.legend()
plt.xlabel("z")
plt.ylabel("r")
plt.title("Comparison between proplyd and isotropic inner winds".format(innertype))
plt.axes().set_aspect("equal")
plt.savefig("shell-test-compare.png")
plt.clf()
