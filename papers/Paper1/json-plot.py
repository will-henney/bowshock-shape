"""
Plot the ouputs from the program
space-param-tab.py
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import json
import argparse
import sys

### Command line section ###

parser = argparse.ArgumentParser(description="Plot output options")

parser.add_argument("-x",type=str,default="i",choices=["i","I","R0","r0","beta","b"],help="x axis quantity")
parser.add_argument("-y",type=str,default="Rc",choices=["Rc","rc","R0","r0"],help="y axis quantity")
#the choices are made in a such way that does not distinguish between uppercases and lowercases, also beta and b
#are equivalent choices

args= parser.parse_args()

inc_choice = (args.x == "i") or (args.x == "I")
r0_choice_x = (args.x == "r0") or (args.x == "R0")
r0_choice_y = (args.y == "r0") or (args.y == "R0")
beta_choice = (args.x == "beta") or (args.x == "b")
rc_choice = (args.y == "Rc") or (args.y == "rc")
print("You chose to plot {} vs {}".format(args.x,args.y))

print(r0_choice_x,"\n",r0_choice_y,"\n",inc_choice,"\n",beta_choice,"\n",rc_choice)
#chosing R0 as x axis and y axis should lead to the identity function
if r0_choice_x and r0_choice_y:
    print("These choices don't lead to useful information")
    sys.exit()

### Open json file according to the x and y axis choices ###

# RcVsR0.json has high resolution in inclinations
# Rc-R0VsBeta.json has high resolution in beta

# The potential plots we can do are the following:

# Requires to open RcVsR0.json
# - Rc vs R0 for a sample of beta values
# - Rc vs inclination for a sample of beta values
# - R0 vs inclination for a sample of beta values

# Requires to open Rc-R0VsBeta.json
# - Rc vs beta for a sample of inclinations
# - R0 vs beta for a sample of inclinations
# - Rc vs R0 for a sample of inclinations


### Plotting section ###

params = {
    "font.family": "serif",
    "text.usetex": True,
    "text.latex.preamble": [r"\usepackage[varg]{txfonts}"],
    "figure.figsize": (5, 5),
}
matplotlib.rcParams.update(params)


if r0_choice_x:
    shelldata = json.load(open("RcVsR0.json"))
    shelldata1 = json.load(open("Rc-R0VsBeta.json"))

    for b, modeldata in shelldata.items():
        inclinations = shelldata[b]["inc"]
        R0 = np.array(shelldata[b]["R0'"])
        Rc = np.array(shelldata[b]["Rc'"])
        plt.plot(R0,Rc/R0,label=r"$\beta={}$".format(b))
    plt.legend(fontsize="small")
    plt.xlabel(r"$R'_0/D'$")
    plt.ylabel(r"$R'_c/R'_0$")
    plt.savefig("Rc-R0-b.pdf")
    plt.clf()

    for i, modeldata in shelldata1.items():
        beta = np.array(shelldata1[i]["beta"])
        R0 = np.array(shelldata1[i]["R0'"])
        Rc = np.array(shelldata1[i]["Rc'"])
        plt.plot(R0,Rc/R0,label=r"$i={}$".format(i))
    plt.legend(fontsize="small")
    plt.xlabel(r"$R'_0/D'$")
    plt.ylabel(r"$R'_c/R'_0$")
    plt.savefig("Rc-R0-i.pdf")

elif inc_choice:
    shelldata = json.load(open("RcVsR0.json"))
    for b, modeldata in shelldata.items():
        inclinations = np.array(shelldata[b]["inc"])
        R0 = np.array(shelldata[b]["R0'"])
        Rc = np.array(shelldata[b]["Rc'"])
        if r0_choice_y:
            Y=R0
        elif rc_choice:
            Y=Rc/R0
        plt.plot(inclinations,Y,label=r"$\beta={}$".format(b))
    plt.legend(fontsize="small")
    plt.xlabel("inclination (deg)")
    if r0_choice_y:
        ylabel=r"$R'_0/D$"
        out = "R0-i.pdf"
    elif rc_choice:
        ylabel=r"$R'_c/D'$"
        out="Rc-i.pdf"
    plt.ylabel(ylabel)
    plt.savefig(out)
    
elif  beta_choice:
    shelldata = json.load(open("Rc-R0VsBeta.json"))
    for i, modeldata in shelldata.items():
        beta = np.array(shelldata[i]["beta"])
        R0 = np.array(shelldata[i]["R0'"])
        Rc = np.array(shelldata[i]["Rc'"])
        if r0_choice_y:
            Y=R0
        else:
            Y=Rc/R0
        plt.plot(beta,Y,label=r"$i={}$".format(i))
    plt.legend(fontsize="small")
    plt.xlabel(r"$\beta$")
    if r0_choice_y:
        ylabel=r"$R'_0/D$"
        out = "R0-b.pdf"
    elif rc_choice:
        ylabel=r"$R'_c/D'$"
        out="Rc-b.pdf"
    plt.ylabel(ylabel)
    plt.savefig(out)


print("Program finished")
