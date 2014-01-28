import numpy as np
import glob
import json
import matplotlib.pyplot as plt
import argparse
"""
This program is designed to capture some important parameters
from other files (R0,Rh,beta,inc) to plot one versus other
"""

parser = argparse.ArgumentParser(description="What vs what I want to plot?")
parser.add_argument("--yaxis",type=str,default="Rh",choices={"Rc","Rh","R90"},help="y axis quantity")
#parser.add_argument("--xaxis",type=str,default="R0",choices={"R0","beta"},help="x axis quantity")
parser.add_argument("--norm",action="store_true",help="normalize y axis with R0")
cmd_args = parser.parse_args()
yaxis = cmd_args.yaxis
norm = cmd_args.norm 

#1
#Find all the json files with the bowshocks info
pattern = ["isotropic-beta-0.???-arcdata.json","proplyd-beta-0.???-arcdata.json"]
inclinations = [0.0,15.0,30.0,45.0,60.0]

#2 Read the json files
for pat in pattern:
    for j in glob.glob(pat):
        f = json.load(open(j))
#3 Extract the relevant information: beta, inclinations, R0, Rh, R90
        beta = f["info"]["beta"]
        for inc in inclinations:
            try:
                Y = f["outer"]["i="+str(inc)][yaxis]
                R0 = f["outer"]["i="+str(inc)]["R0"]
            except KeyError:
                break
                #maybe for high beta values the desired inclination does not exist
            
#4 Plot everything vs everything, but being able to choose only one plot per run 
# I want Rc,Rh,R90 over D or over R0 Vs R0/D or beta
            
            if norm:
                plt.plot(R0,Y/R0,"r.")
            else:
                plt.plot(R0,Y,"g.")

plt.grid()
plt.xlabel("R0'/D'")
plt.ylabel("y axis")
plt.title("Title")
plt.savefig("something_vs_r0.pdf")
