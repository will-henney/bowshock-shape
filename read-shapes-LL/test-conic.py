import numpy as np
import matplotlib.pyplot as plt
import conic_utils

def x90(R0,th):
    """
    R90 value in function of R0, it is kind complicated
    """
    B = np.tan(np.radians(th))
    s = np.sign(B)
    if B ==0:
        x = np.sqrt(2*R0)
    else:
        x = np.sqrt((s/B**2)*( 1-B**4*( -R0+s/B**2 )**2 ) )   
    return x

R0 = 0.3
x = np.linspace(-2.0, 2.0, 1000)
for th in [-75.0, -60.0, -45.0, -30.0, -15.0, 0.0, 15.0, 30.0, 45.0, 60.0]:
    plt.plot(x, conic_utils.yconic_th(x, th))
    plt.plot([x90(R0,th)],[R0],"o")
plt.axis('equal')
plt.xlim(-2, 2)
plt.ylim(0, 2)
plt.savefig("test/test-conic.pdf")
