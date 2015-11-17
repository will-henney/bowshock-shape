import numpy as np
import json
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import argparse

"""
Attempt to fit log beta vs log [Rc,R0,R90]
"""

shelldata = json.load(open("Rc-R0VsBeta.json"))

parser = argparse.ArgumentParser(description="choose y-axis to plot")
parser.add_argument("--axis",type=str,default="Rc",choices=["Rc","R0","R90"],help="y-axis choice")
parser.add_argument("--norm",action="store_true",help="Normalize with R0")

cmd_args = parser.parse_args()
b = np.logspace(-3,-1)
y_arg = {"R0":"R0'","Rc": "Rc'_f","R90":"R90'"}
y_lab = {"R0":"R'_0","Rc": "R'_c","R90":"R'_{90}"}
y_a = {"R0":0.5*np.log10(b)-np.log10(1+np.sqrt(b)),"Rc":np.log10(1.5)+0.5*np.log10(b)-np.log10(1-b),
       "R90":0.5*(np.log10(2.4)+np.log10(b)-3*np.log(1-0.8*b))}
y_n = {"R0":0*b,"Rc":np.log10(1.5) - np.log10(1.-np.sqrt(b)) ,"R90":0.5*(np.log10(2.4)+2*np.log10(1+np.sqrt(b))-3*np.log(1-0.8*b))}
for i,modeldata in shelldata.items():
    beta = np.array(shelldata[i]["beta"])
    y = np.array(shelldata[i][y_arg[cmd_args.axis]])
    if cmd_args.norm:
        R0 = np.array(shelldata[i]["R0'"])
        Y = y/R0
    else:
        Y = y
    plt.plot(np.log10(beta),np.log10(Y),":",label=i)

if cmd_args.norm:
    ylabel = r"$\log\left({}/R'_0\right)$".format(y_lab[cmd_args.axis])
    fig_name = "log-beta-vs-log-{}-norm.pdf".format(y_lab[cmd_args.axis])
    y_analytic = y_n[cmd_args.axis]
else:
    ylabel = r"$\log{}$".format(y_lab[cmd_args.axis])
    fig_name = "log-beta-vs-log-{}.pdf".format(y_lab[cmd_args.axis])
    y_analytic = y_a[cmd_args.axis]
plt.plot(np.log10(b),y_analytic,"k--",lw=2,alpha=0.5,label="analytic i=0")
plt.legend(loc="best",fontsize="small")
plt.xlabel(r"$\log\beta$")
plt.ylabel(ylabel)
plt.savefig(fig_name)
