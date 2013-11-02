import numpy as np
import matplotlib.pyplot as plt
import conic_utils

x = np.linspace(-2.0, 2.0, 1000)
for th in [-75.0, -60.0, -45.0, -30.0, -15.0, 0.0, 15.0, 30.0, 45.0, 60.0]:
    plt.plot(x, conic_utils.yconic_th(x, th))
plt.axis('equal')
plt.xlim(-2, 2)
plt.ylim(0, 2)
plt.savefig("test/test-conic.pdf")
