from __future__ import print_function
import numpy as np
import sys
sys.path.insert(0,"../../conic-projection")
import conproj_utils as conproj
import matplotlib.pyplot as plt

"""
This scipt has been created for making some figures in section
4 
"""

############## Variying  \theta_c with unitary radius of curvature ####################

#45: circle
#15: ellipse
# 0: parabola
#-15: hyperbola
#-45:hyperbola

for t in [45,15,-15,-45]:
    c  = conproj.Conic(A=1.0,th_conic = t)
    trange = c.make_t_array()
    plt.plot(c.x(trange),c.y(trange),label =r"$\theta_c={}^\circ$".format(t))

plt.legend(loc="best")
plt.xlim(-3,1.1)
plt.gca().set_aspect("equal",adjustable="box") #Trick to set equal axes found in stackoverflow.com
plt.ylim(-2.1,2.1)
plt.xlabel(r"$x/R_0$")
plt.ylabel(r"$y/R_0$")
plt.savefig("conic1.pdf")

##########  Show same conic at different inclinations (multiple examples) #############


# i=0
# i=15
# i=30
# i=45
