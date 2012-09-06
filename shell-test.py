import numpy as np 
from equation6 import Shell
import matplotlib.pyplot as plt

N = 20                         # number of angles
theta = np.linspace(0.0, np.pi, N)

# Test isotropic shells
for beta in 1.0, 0.2, 0.05, 0.02, 0.005:
    shell = Shell(beta=beta)
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
plt.axes().set_aspect("equal")
plt.savefig("shell-test-isotropic.png")
