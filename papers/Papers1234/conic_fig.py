"""
This scipt has been created for making some figures in section 4 
"""
from __future__ import print_function
import numpy as np
import sys
sys.path.insert(0,"../../conic-projection")
import conproj_utils as conproj
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


parser = argparse.ArgumentParser(description="Choose the output figure")
parser.add_argument("--fig",type=int,default = 0,choices=[0,1],help="Choose your figure output")
args = parser.parse_args()
flag= args.fig <= 0 
############## Variying  \theta_c with unitary radius of curvature ####################

# Figure style
sns.set_style('ticks')
#sns.set_context('talk')
sns.set_color_codes('dark')


#45: circle
#30-15: ellipse
# 0: parabola
#-15- -45:hyperbola

thc_list = [60, 45, 30, 0.0, -30,-45]
conic_types = ['Oblate ellipse', 'Circle', 'Prolate ellipse',
               'Parabola', 'Hyperbola', 'Hyperbola']

if flag:
    fig, ax = plt.subplots(1, 1)
    ax.axhline(xmin=0.05,xmax=0.95,color="black",lw=1.5,alpha=0.7)
    ax.axvline(ymin=0.05,ymax=0.95,color="black",lw=1.5,alpha=0.7)
    
    for t, s in reversed(list(zip(thc_list, conic_types))):
        c  = conproj.Conic(A=2.0,th_conic = t)
        trange = c.make_t_array()
        ax.plot(c.x(trange), c.y(trange), lw=3, alpha=0.8,
                label =r"{}: $\theta_q={:.0f}^\circ$".format(s, t))

    # Star at origin
    ax.plot([0],[0], "w*", ms=17, alpha=0.9)
    ax.plot([0],[0], "b*", ms=15)
    ax.plot([0],[0], "w*", ms=10, alpha=0.95)
    # Circle at apex
    ax.plot([1],[0], "wo", ms=7, alpha=0.9)
    ax.plot([1],[0], "bo", ms=6)
    ax.plot([1],[0], "wo", ms=4, alpha=0.95)

    legend = ax.legend(loc="lower left", frameon=True, fontsize='small')
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(0.95)
    ax.set_xlim(-5,2.1)
    ax.set_aspect("equal", adjustable="box") # Set equal axis scales in x, y
    ax.set_ylim(-4.1,4.1)
    ax.set_xlabel(r"$x/R_0$")
    ax.set_ylabel(r"$y/R_0$")
    fig.set_size_inches(5, 5)
    fig.tight_layout()
    sns.despine()
    fig.savefig("conic1.pdf")

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
# plot y'/R'_0 vs x'/R'_0
        ax.plot(ch.xt(i,t)/(ch.g(i)*np.cos(np.radians(i))),ch.yt(i,t)/(ch.g(i)*np.cos(np.radians(i))))
        Ap = ch.Aprime(i)
        ax.plot(Ap*np.cos(t)-Ap+1,Ap*np.sin(t),label=r"A'={:.2f}".format(Ap))
        ax.plot([0],[0],"b*")
        ax.legend(loc="best",fontsize="small")
        plt.axhline(xmin=0.05,xmax=0.95,color="black",lw=1.5,alpha=0.7)
        plt.axvline(ymin=0.05,ymax=0.95,color="black",lw=1.5,alpha=0.7)
        ax.set_xlabel(r"x'/R'_0")
        ax.set_ylabel(r"y'/R'_0")
        ax.set_xlim(-3,1.5)
        ax.set_ylim(-2.1,2.1)
        ax.set_title("i={}".format(i))
    f.set_size_inches(10,10)
    f.savefig("conic2.pdf")
