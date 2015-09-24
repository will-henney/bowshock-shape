"""
Create a json dictionary for one of the following figures
0. Rc', R0' vs beta for i=0,15,30,45,60 degrees
1. Json dictionary for Rc' vs R0' for beta = 0.001,0.01,0.05,0.08,0.1 if possible
"""

import numpy as np
import sys
sys.path.insert(0,"../..")
from equation6 import Shell
import toolbox as tool
import matplotlib.pyplot as plt
import argparse
import json

parser = argparse.ArgumentParser(description="Terminal options for script")
parser.add_argument("--test",action="store_true",help="Run test for Rc accuracy")
parser.add_argument("--fig",type=int,default=0,choices=[0,1],help="choose figure output")
args= parser.parse_args()

theta = np.linspace(0,np.pi,500) #Looks like 500 points are a good enough resolution to
                                 #obtain good results
#Test Rc function accuracy

if args.test:
    shell = Shell(beta=0.01,innertype="proplyd")
    R = shell.radius(theta)
    R[R<=0]=np.nan
    Rc = tool.Rc(R,theta)
    Rc_analytic = 1.5*R[0]/(1-np.sqrt(0.01))
    print(Rc,Rc_analytic,Rc-Rc_analytic)
    print("Test finished")
    sys.exit() #After running the test exit python

if args.fig==0:
    inc=[0,15,30,45,60]
    beta = np.logspace(-3,-1,30)
    shelldata = {}
    for i in inc:
        print("Case i={}".format(i))
        R0p = np.zeros_like(beta)*np.nan
        Rcp = np.zeros_like(beta)*np.nan

        for j,b in enumerate(beta):
            if i==0: #Avoid recalculate R for each iteration of inclination
                shell = Shell(beta=b,innertype="proplyd")
                R = shell.radius(theta)
                R[R<=0]=np.nan
            xp,yp = tool.projection(R,theta,i)
            Rp = tool.R_prime(xp,yp)
            tp = tool.theta_prime(xp,yp)
            try:
                R0p[j] = Rp[0]
                Rcp[j] = tool.Rc(Rp,tp)
            except IndexError:
                break

        shelldata[i]={
            "beta": beta[np.isfinite(R0p)].astype(float).tolist(),
            "R0'":  R0p[np.isfinite(R0p)].astype(float).tolist(),
            "Rc'":  Rcp[np.isfinite(R0p)].astype(float).tolist()
        }
    with open("Rc-R0VsBeta.json","w") as f:
        json.dump(shelldata,f,indent=2)
    
if args.fig ==1:
    beta=[0.001,0.01,0.05,0.08,0.1]
    inc = np.linspace(0,60)
    shelldata = {}
    for b in beta:
        print("Creating shell with beta={}".format(b))
        shell = Shell(beta=b,innertype="proplyd")
        R= shell.radius(theta)
        R[R<=0] = np.nan
        R0p = np.zeros_like(inc)*np.nan
        Rcp = np.zeros_like(inc)*np.nan

        for j,i in enumerate(inc):
            xp,yp = tool.projection(R,theta,i)
            Rp = tool.R_prime(xp,yp)
            tp = tool.theta_prime(xp,yp)
            try:
                R0p[j]=Rp[0]
                Rcp[j]=tool.Rc(Rp,tp)
            except IndexError:
                print("Maximum inclination {:.2f}".format(i))
                break
        
        shelldata[b] = {
            "inc": inc[np.isfinite(R0p)].astype(float).tolist(),
            "R0'": R0p[np.isfinite(R0p)].astype(float).tolist(),
            "Rc'": Rcp[np.isfinite(R0p)].astype(float).tolist()}

    with open("RcVsR0.json","w") as f:
        json.dump(shelldata,f,indent=2)
