
import numpy as np
import json
from scipy.interpolate import interp1d

#####################################################
# My first program written and tangled with org mode#
#####################################################

# Selecting the wanted *-arcdata.json files

path_to_files = "../LL-shapes-2017/"
files_pattern = "-arcdata.json"
ID = [
     "w073-227","109-246","w000-400",
     "w005-514","w012-407","w030-524",
     "042-628","LL1","065-502",
     "w069-601","4285-458","LL3",
     "LL2","LL4","4468-605",
     "116-3101","w266-558","308-3036",
     "LL6"
]
# Extract desired info: R, theta, R0, Rc for inner and outer shell
R_out = {}
theta_out = {}
Rc_out = {}
R0_out = {}
R_in = {}
theta_in = {}
Rc_in = {}
R0_in = {}
#Set R90 dictionaries 
R90_out = {}
Rm90_out = {}
R90_in = {}
Rm90_in = {}

for id in ID:
    file = path_to_files + id + files_pattern
    data = json.load(open(file))
    R_out[id] = data["outer"]["R"]    
    theta_out[id] = data["outer"]["theta"]
    Rc_out[id] = data["outer"]["Rc"]
    R0_out[id] = data["outer"]["R0"]
    try:
        R_in[id] = data["inner"]["R"]
    except KeyError:
        R_in[id] = None
    try:   # Looks like not all objects have inner shell data
        theta_in[id] = data["inner"]["theta"]
    except KeyError:
        theta_in[id] = None
    try:
        Rc_in[id] = data["inner"]["Rc"]
    except KeyError:
        Rc_in[id] = None
    try:
        R0_in[id] = data["inner"]["R0"]    
    except KeyError:
        R0_in[id] = None

    # Interpolate R90 and/or Rm90 (the opposite R90) from data
    # Outer shell
    Ro_int = interp1d(np.array(theta_out[id]),np.array(R_out[id]))
    if max(theta_out[id]) >= 90.0:
        R90_out[id] = float(Ro_int(90.0))
    else:
        R90_out[id] = None
    if min(theta_out[id]) <= -90.0:
        Rm90_out[id] = float(Ro_int(-90.0))
    else:
        Rm90_out[id] = None 
    # Inner shell
    if theta_in[id] is not None and R_in[id] is not None:
        Ri_int = interp1d(np.array(theta_in[id]),np.array(R_in[id]))
        if max(theta_in[id]) > 90.0:
            R90_in[id] = float(Ri_int(90.0))
        else:
            R90_in[id] = None
        if min(theta_in[id]) < -90.0:
            Rm90_in[id] = float(Ri_int(-90.0))
        else:
            Rm90_in[id] = None
    else:
        R90_in[id] = None
        Rm90_in[id] = None
# Put all the data in a single dictionary
output = {
         "outer":{
                 "R0":R0_out,
                 "Rc":Rc_out,
                 "R90":R90_out,
                 "Rm90":Rm90_out
    },
         "inner":{
                 "R0":R0_in,
                 "Rc":Rc_in,
                 "R90":R90_in,
                 "Rm90":Rm90_in
    },
         "help":{
                "outer":"Means outer shell",
                "inner":"Means inner shell",
                "R0":"Shell radius along symmetry axis",
                "Rc":"Radius of curvature at theta=0",
                "R90":"Shell radius perpendicular to R0",
                "Rm90":"Same as R90 but in the opposite direction"
    }
}
# Write the output in a json file
with open("radii-set.json","w") as f:
    json.dump(output,f,indent=4)
