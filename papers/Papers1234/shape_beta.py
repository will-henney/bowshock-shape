import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,"../..")
from equation6 import Shell
import argparse
import seaborn as sns
###################################
"""                              
The objective of this program is 
#to make a figure that illustrates
how the shell shape changes with
$\beta$ and if possible compare it
with the result of other authors
"""
###################################

############ Figure 1 #############
# Shell shape vs $\beta$ ##########

beta = [0.01,0.1,0.99]
theta = np.linspace(0,0.99*np.pi)

parser =argparse.ArgumentParser(description="choose figure output")
parser.add_argument("--fig",type=int,default=0,choices=[0,1],help="figure output")
args = parser.parse_args()
flag = args.fig <= 0
inner_list = ["isotropic","proplyd"]
xi_list = [1.0,0.8,0.2]
line_dict = {"isotropic":"--","proplyd":"-"}
line_list = ["--","-",":",]
color_list = ["r","g","b","m","c","k"]
#label = {"isotropic":None,"proplyd":r"$\beta={}$".format(b)}
sns.set_style("whitegrid")
if flag:
    for xi,line in zip(xi_list,line_list):
        for b,col in zip(beta,color_list):
            shell = Shell(beta=b,innertype="anisotropic",xi=xi)
            R = shell.radius(theta)
            R[R<=0] = np.nan
            if xi==0.8:
                label = r"$\beta={}$".format(b)
            else:
                label = None
            plt.plot(R*np.cos(theta),R*np.sin(theta),color=col,linestyle=line,
                     label=label)
    fontsize = 15
    ticksize = 14
    plt.legend(fontsize="small")
    plt.xlabel(r"$z/D$",fontsize=fontsize)
    plt.ylabel(r"$r/D$",fontsize=fontsize)
    plt.tick_params(axis='both', which='major', labelsize=ticksize)
    plt.gca().set_aspect("equal",adjustable="box")
    fig = plt.gcf()
    fig.set_size_inches(7, 4.5)
    plt.xlim(-0.4,1)
    plt.ylim(0,1)
    plt.tight_layout()
    fig.savefig("figs/r-beta.pdf")

####### Figure 2 #################

# Show the characteristic Radii for a generic bowshock


else:
    b = 0.1
    t = np.linspace(-np.pi,np.pi)
    BS = Shell(beta=b,innertype="proplyd")
    A = 1.5/(1-np.sqrt(b))
    R = BS.radius(theta)
    R[R<=0]=np.nan
    plt.plot(R*np.cos(theta),R*np.sin(theta),"k-",R*np.cos(theta),-R*np.sin(theta),"k-",lw=3)
    plt.plot(A*R[0]*np.cos(t)-R[0]*(A-1),A*R[0]*np.sin(t))
    plt.plot([0],[0],"r*")
    plt.plot([-R[0]*(A-1)],[0],"k.")
    plt.grid()
    plt.gca().set_aspect("equal",adjustable="box")
    plt.xlabel("z/D")
    plt.ylabel("r/D")
    plt.xlim(-R[0]*(A+1)-0.1,R[0]+0.6)
    plt.ylim(-0.3,1.5)
    plt.savefig("ch-radii.pdf") 
    














