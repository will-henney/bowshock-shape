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
parser.add_argument("-y",type=str,default="Rc",choices=["Rc","rc","R0","r0","r90","R90"],help="y axis quantity")
#the choices are made in a such way that does not distinguish between uppercases and lowercases, also beta and b
#are equivalent choices

args= parser.parse_args()

inc_choice = (args.x == "i") or (args.x == "I")
r0_choice_x = (args.x == "r0") or (args.x == "R0")
r0_choice_y = (args.y == "r0") or (args.y == "R0")
beta_choice = (args.x == "beta") or (args.x == "b")
rc_choice = (args.y == "Rc") or (args.y == "rc")
r90_choice = (args.y == "R90") or (args.y == "r90")
print("You chose to plot {} vs {}".format(args.y,args.x))

#print(r0_choice_x,"\n",r0_choice_y,"\n",inc_choice,"\n",beta_choice,"\n",rc_choice)
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
    color = "rgbmyck"
    if rc_choice:
   #     Y = Rc
        ylab =r"$R'_c/R'_0$"
        out = "Rc-R0-b.pdf"
        out1 = "Rc-R0-i.pdf"
        
        
    elif r90_choice:
    #    Y = R90
        ylab = r"$R'_{90}/R'_0$"
        out = "R90-R0-b.pdf"
        out1 = "R90-R0-i.pdf"
    iteration=0
    for b, modeldata in shelldata.items():
        inclinations = modeldata["inc"]
        R0 = np.array(modeldata["R0'"])
        Rc = np.array(modeldata["Rc'_f"])
        R90 = np.array(modeldata["R90'"])
        every15 = np.zeros(R0.shape,dtype=bool)
        for inc15 in 0.0,15.0,30.0,45.0,60.0,75.0:
            iclosest = np.argmin(np.abs(np.array(inclinations)-inc15))
            every15[iclosest]=True
        if rc_choice:
            Y=Rc/R0
            #dY = np.array(shelldata[b]["dRc'"])
            plt.plot(R0,Y,lw=3,c=color[iteration],label=r"$\beta={}$".format(b))
            plt.plot(R0[every15],Y[every15],"d",c=color[iteration])
        else:
            Y=R90/R0
            plt.plot(R0,Y,c=color[iteration],lw=3,label=r"$\beta={}$".format(b))
            plt.plot(R0[every15],Y[every15],"d",c=color[iteration])
        iteration+=1
        if iteration>6:
            iteration-=7
    plt.legend(fontsize="small")
    plt.xlabel(r"$R'_0/D'$")
    plt.ylabel(ylab)
    plt.xlim(-1e-4,0.5+1e-4)
    plt.ylim(-1e-4,5+1e-4)
    plt.savefig(out)
    plt.clf()
    iteration=0
    for i, modeldata in shelldata1.items():
        beta = np.array(shelldata1[i]["beta"])
        R0 = np.array(shelldata1[i]["R0'"])
        Rc = np.array(shelldata1[i]["Rc'_f"])
        R90 = np.array(shelldata1[i]["R90'"])
        bpoint = np.zeros(R0.shape,dtype=bool)
        for bskip in 0.002,0.005,0.01,0.05,0.1:
            bclosest = np.argmin(np.abs(beta-bskip))
            bpoint[bclosest] = True

        if rc_choice:
            Y=Rc/R0
            #dY = np.array(shelldata1[i]["dRc'"])
            plt.plot(R0,Y,c=color[iteration],lw=3,label=r"$i={}$".format(i))
            plt.plot(R0[bpoint],Y[bpoint],"o",c=color[iteration])
        else:
            Y=R90/R0
            plt.plot(R0,Y,c=color[iteration],lw=3,label=r"$i={}$".format(i))
            plt.plot(R0[bpoint],Y[bpoint],"o",c=color[iteration])
        iteration+=1
        if iteration>6:
            iteration-=7
    plt.legend(fontsize="small")
    plt.xlabel(r"$R'_0/D'$")
    plt.ylabel(ylab)
    plt.xlim(-1e-4,0.5+1e-4)
    plt.ylim(-1e-4,5+1e-4)
    plt.savefig(out1)

elif inc_choice:
    shelldata = json.load(open("RcVsR0.json"))
    for b, modeldata in shelldata.items():
        inclinations = np.array(shelldata[b]["inc"])
        R0 = np.array(shelldata[b]["R0'"])
        Rc = np.array(shelldata[b]["Rc'_f"])
        R90 = np.array(shelldata[b]["R90'"])
        if r0_choice_y:
            Y=R0
            ylabel=r"$R'_0/D$"
            out = "R0-i.pdf"
            dy = None
        elif rc_choice:
            Y=Rc
            ylabel=r"$R'_c/D'$"
            out="Rc-i.pdf"
           # dy = np.array(shelldata[b]["dRc'"])
        elif r90_choice:
            Y = R90
            ylabel=r"$R'_{90}/D'$"
            out="R90-i.pdf"
            #dy = None
        plt.plot(inclinations,Y,label=r"$\beta={}$".format(b))
    plt.legend(fontsize="small")
    plt.xlabel("inclination (deg)")
    plt.ylabel(ylabel)
    plt.savefig(out)
    
elif  beta_choice:
    shelldata = json.load(open("Rc-R0VsBeta.json"))
    for i, modeldata in shelldata.items():
        beta = np.array(shelldata[i]["beta"])
        R0 = np.array(shelldata[i]["R0'"])
        Rc = np.array(shelldata[i]["Rc'_f"])
        R90 = np.array(shelldata[i]["R90'"])
        if r0_choice_y:
            Y=R0
            ylabel=r"$R'_0/D'$"
            out = "R0-b.pdf"
            y_analytic = np.sqrt(beta)/(1+np.sqrt(beta))
            #dy=None
        elif rc_choice:
            Y=Rc
            ylabel=r"$R'_c/D'$"
            out = "Rc-b.pdf"
            y_analytic = R0/(1-np.sqrt(beta))
            #dy = np.array(shelldata[i]["dRc'"])
        elif r90_choice:
            Y = R90
            ylabel=r"$R'_{90}/D'$"
            out = "R90-b.pdf"
            y_analytic = np.sqrt(2.4*beta/(1-0.8*beta)**3)
            #dy=None
        plt.plot(beta,Y,label=r"$i={}$".format(i))
    plt.plot(beta,y_analytic,"k--",label = "analytic i=0")
    plt.legend(fontsize="small")
    plt.xlabel(r"$\beta$")
    plt.ylabel(ylabel)
    plt.savefig(out)


print("Program finished")
