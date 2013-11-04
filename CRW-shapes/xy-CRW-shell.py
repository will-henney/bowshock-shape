from equation6 import Shell 
import numpy as np
import astropy
from astropy import coordinates as coord
from region_utils import region_point_to_string
import argparse
import json
from misc_utils import run_info
"""
Goal: The idea of this program is create a CRW bowshock and write an output *-forma.reg file to use it later with Will's former programs
Steps: 
1.- Create the bowshock arc with the equation6 library for a fixed beta parameter for isotropic and/or proplyd case (perform to any beta value for both 
cases).
2.- Assign a fictional RA-Dec proplyd position (maybe the same position of a known proplyd)
3.- Use the astropy packages as well other libraries to write the coordinates of the arc in a ds9 region format  
"""

#1
beta = 0.01
theta = np.linspace(0,0.5*np.pi,500)
shell = Shell(beta,"isotropic")
R = shell.radius(theta)

arcdata = {
    "info":{
        "Description": "CRW theoretical stationary bowshocks",
        "history":["Initially created by "+ run_info()],
        "beta": beta,
        "i": 0.0,
        "inner wind type": "isotropic"
        },
    "outer":{
        "x": list(R*np.cos(theta)),
        "y": list(R*np.sin(theta))
        },
    "help":{
        "beta":"Two-winds momentum ratio",
        "i": "Rotation angle between proplyd reference frame and observer reference frame",
        "inner wind  type": "Angular variation of density's inner wind",
        "x": "(list) x-coordinates of the outer shell",
        "y": "(list) y-coordinates of the outer shell",
        }

}

jsonfile= "teo-shell-beta-{}.json".format(beta)
with open(jsonfile,"w") as f:
    json.dump(arcdata,f,indent=4)


