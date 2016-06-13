import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm
import argparse
import seaborn as sns
import sys
sys.path.insert(0,"../conic-projection")
from conproj_utils import Conic
import glob
import json

# Contrast  conic curves against observations


def A(b, xi):
    """
    Returns the radius of curvature normalized with R0

    Corrected version that also depends on anisotropy index k

    xi = 2/(k+2) => k = 2 (1/xi - 1)
    """
    k = 2*(1./xi - 1.)
    sb = np.sqrt(b)
    c = (1 - b - 9.0*k/4.0)/30.0
    alpha = (1.0 + 2*sb)/6.0 + c/(1 + sb)
    return 1./(1.0-2*alpha)

def theta_c(beta,xi=1.0):
    """
    theta_c defines the excentricity of a given conic
    """
    Acurv = A(beta,xi)
    arg = 2*Acurv - 3*xi*(1.0 + np.sqrt(beta))**2/(1.0 - xi*beta)**2/(1 + 0.2*xi*beta)
    return np.sign(arg)*np.arctan(np.sqrt(np.abs(arg)))

def q(b):
    return np.sqrt(b)/(1.+np.sqrt(b))
#Step 1 Define command variables

parser = argparse.ArgumentParser(description="Command line options")
parser.add_argument("--xi",type=float,default=1.0,help="inner wind parameter")

args = parser.parse_args()
Xi = args.xi

#Step 2: Initialize arrays


beta = [5e-4, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1]
sns.set_style("ticks")
colors = sns.color_palette('Blues', len(beta))

inc = np.linspace(0,75,100)
#print(inc)

#step 3: Beta loop  for plotting
for b, c in zip(beta, colors):
    conic = Conic(A=A(b, Xi),th_conic=np.degrees(theta_c(b,Xi)))
    print(r"$\beta={}$".format(b))
    q_prime = q(b)*conic.g(inc)
    A_prime = conic.Aprime(inc)
    plt.plot(q_prime,A_prime, linestyle="-", linewidth=2.0, color=c,
	label=r"$\beta={}${}$\theta_c={:.0f}^\circ$".format(b,"\n",np.degrees(theta_c(b,Xi))))
    every15 = np.zeros(q_prime.shape, dtype=bool)
    for thisinc in 0.0, 15.0, 30.0, 45.0, 60.0, 75.0:
        iclosest = np.argmin(np.abs(inc - thisinc))
        every15[iclosest] = True
        plt.plot(q_prime[every15], A_prime[every15], '.', color=c)
#step 4: Add the observational points


#collection of hex colors
dark_blue = "#1e25b6"
pearl_turquoise ="#32c6a6"
mexican_pink = "#e4007c"
crimson = "#dc143c"
leaf_green = "#15ae26"
brown = "#b6451e"
gray = "#515952"
guinda = "#aa1c47"
gold = "#FFD700"
orange = "#E08000"
#Create a dictionary with hex colors for the objects
colordict = {"LV2":dark_blue,"LV2b":pearl_turquoise,"LV3":mexican_pink,"LV4":crimson,"LV5":brown,
             "168-328":leaf_green,"169-338":gray,"177-341":guinda,"180-331":orange}

m_savefiles = glob.glob("LV-bowshocks-xyfancy-positionswill-*.save")
dict_xtext = {"LV2":10,"LV2b":-10,"LV3":10,"LV4":10,"LV5":10,"168-328":-10,"169-338":-10,"177-341":10,"180-331":-20}
dict_ytext = {"LV2":10,"LV2b":10,"LV3":-10,"LV4":10,"LV5":10,"168-328":10,"169-338":10,"177-341":-10,"180-331":-20}
for savefile in m_savefiles:
    data = json.load(open(savefile))
    combined_file = savefile.replace('positionswill', 'variations')
    vardata = json.load(open(combined_file))
    plt.plot(data["R0"],data["Rc"]/data["R0"],
             # color=colordict[data["proplyd"]],
             color='k',
             marker="o")
    plt.annotate(data["proplyd"], xy=(data["R0"],data["Rc"]/data["R0"]),
                 xytext=(dict_xtext[data["proplyd"]], dict_ytext[data["proplyd"]]),
                 textcoords="offset points", fontsize="xx-small",
                 bbox=dict(boxstyle='round,pad=0.5',
                           fc=colordict[data["proplyd"]],
                           alpha=0.5))
    # Plot the variations of the fits with points removed
    R0 = data["R0"]
    A = data["Rc"]/data["R0"]
    var_R0 = vardata["R0"]
    var_A = np.array(vardata["Rc"])/np.array(vardata["R0"])
    for vR0, vA in zip(var_R0, var_A):
        # Scale gives fractional deviation from typical value
        scale = np.hypot((vR0 - R0)/0.25, (vA - A)/1.5)
        alpha = 1./(1 + 20.0*scale)
        plt.plot([R0, vR0], [A, vA], '-',
                  lw=2, alpha=alpha, color=colordict[data["proplyd"]])
epsilon=1e-4
plt.grid()
legend = plt.legend(loc="upper left",fontsize="small", ncol=2,
                    fancybox=True, frameon=True)
legend.get_frame().set_facecolor('white')
plt.xlim(epsilon,0.4+epsilon)
plt.ylim(epsilon,4+epsilon)
plt.xlabel(r"$R'_0/D'$")
plt.ylabel(r"$R'_c/R'_0$")
plt.title(r"$\xi={:.1f}$".format(Xi))
plt.savefig("conic_xi-{:02.0f}.pdf".format(10*Xi)) #avoids writing a period in the file name
