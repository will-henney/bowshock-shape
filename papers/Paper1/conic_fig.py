from __future__ import print_function
import numpy as np
import sys
sys.path.insert(0,"../../conic-projection")
import conproj_utils as conproj
import matplotlib.pyplot as plt
import argparse
"""
This scipt has been created for making some figures in section
4 
"""
parser = argparse.ArgumentParser(description="Choose the output figure")
parser.add_argument("--fig",type=int,default = 0,choices=[0,1],help="Choose your figure output")
args = parser.parse_args()
flag= args.fig <= 0 
############## Variying  \theta_c with unitary radius of curvature ####################

#45: circle
#30-15: ellipse
# 0: parabola
#-15- -45:hyperbola
if flag:
    for t in [60, 45, 30, 0.01, -30,-45]:
        c  = conproj.Conic(A=2.0,th_conic = t)
        trange = c.make_t_array()
        plt.plot(c.x(trange),c.y(trange),label =r"$\theta_c={:.0f}^\circ$".format(t))

    plt.legend(loc="center left", fontsize='small')
    plt.xlim(-5,2.1)
    plt.gca().set_aspect("equal",adjustable="box") #Trick to set equal axes found in stackoverflow.com
    plt.ylim(-4.1,4.1)
    plt.xlabel(r"$x/R_0$")
    plt.ylabel(r"$y/R_0$")
    plt.savefig("conic1.pdf")

##########  Show same conic at different inclinations (hyperbola example) #############


# i=0
# i=15
# i=30
# i=45
else:
    f = plt.figure()
    inc = [0,15,30,45]
    ch = conproj.Conic(A=1.0,th_conic = -15)
    t = ch.make_t_array()
    for n,i in enumerate(inc):
        ax = f.add_subplot(2,2,n+1,aspect="equal")
        ax.plot(ch.xt(i,t),ch.yt(i,t))#,label=r"i={}^\circ".format(i))
        Ap = ch.Aprime(i)
        ax.plot(Ap*np.cos(t),Ap*np.sin(t))#,label=r"A'={}".format(Ap))
#        ax.legend(loc="best")
        ax.set_xlabel(r"x/R_0")
        ax.set_ylabel(r"y/R_0")
        ax.set_xlim(-3,1.1)
        ax.set_ylim(-2.1,2.1)
    f.set_size_inches(10,10)
    f.savefig("conic2.pdf")
