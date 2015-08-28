import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,"../..")
from equation6 import Shell
import argparse

###################################
"""                              
The objective of this program is 
#to make a figure that illustrates
how the shaell shape changes with
$\beta$ and if possible compare it
with the result of other authors
"""
###################################

############ Figure 1 #############
# Shell shape vs $\beta$ ##########

beta = [0.001,0.005,0.01,0.05,0.1,0.5,0.99]
theta = np.linspace(0,0.99*np.pi)

parser =argparse.ArgumentParser(description="choose figure output")
parser.add_argument("--fig",type=int,default=0,choices=[0,1],help="figure output")
args = parser.parse_args()
flag = args.fig <= 0

if flag:
    for b in beta:
        shell = Shell(beta=b,innertype="proplyd")
        R = shell.radius(theta)
        R[R<=0] = np.nan
        plt.plot(R*np.cos(theta),R*np.sin(theta),label=r"$\beta={}$".format(b))

    plt.legend(fontsize="small")
    plt.xlabel(r"x/D")
    plt.ylabel(r"y/D")
    plt.gca().set_aspect("equal",adjustable="box")
    plt.xlim(-0.4,1)
    plt.ylim(-1,1)
    plt.savefig("r-beta.pdf")

####### Figure 2 #################

# Show the characteristic Radii for a generic bowshock


else:
    b = 0.1
    t = np.linspace(-np.pi,np.pi)
    BS = Shell(beta=b,innertype="proplyd")
    A = 1./(1-np.sqrt(b))
    R = BS.radius(theta)
    R[R<=0]=np.nan
    plt.plot(R*np.cos(theta),R*np.sin(theta),"k-",R*np.cos(theta),-R*np.sin(theta),"k-")
    plt.plot(A*R[0]*np.cos(t)-R[0]*(A-1),A*R[0]*np.sin(t))
    plt.plot([0],[0],"r*")
    plt.plot([-R[0]*(A-1)],[0],"k.")
    plt.grid()
    plt.gca().set_aspect("equal",adjustable="box")
    plt.xlim(-R[0]*(A+1)-0.1,R[0]+0.1)
    plt.ylim(-0.7,0.7)
    plt.savefig("ch-radii.pdf") 
    














