import sys
sys.path.append("projected")
import bowfuncs as bow
import numpy as np
import matplotlib.pyplot as plt

#print bowfuncs.xt(1.5, 0.1, 0.0)

#steps:
# 0 Array for theta, beta and i
# 1 Create bowshock for different i and beta
# 2 select R0 and R90 for each i
# 3 Make the plot

N = 100
beta = np.array([0.15, 0.05, 0.02, 0.005,0.2,0.3])
i = np.linspace(0,30*np.pi/180.,5)

for b in beta:
    theta = np.linspace(0.01,bow.theta_lim(b),100)
    for j in i:
        xth,yth = np.array([bow.xt(t,b,j) for t in theta ]), np.array([bow.yt(t,b,j) for t in theta ])
        plt.plot(xth,yth,label='i={}'.format(j))
    plt.axis([-1.0,2.0,-0.05,2.0])
    plt.legend()
    plt.xlabel('z')
    plt.ylabel('r')
    plt.title('Bowshock shapes for isotropic wind and beta = {}'.format(b))
    plt.savefig('isotropic-will-wind-plane-sky-beta-{}.png'.format(b))
    plt.clf()
