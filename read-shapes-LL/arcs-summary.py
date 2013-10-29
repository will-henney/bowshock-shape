import glob
import json

pattern = "j8oc??010_drz/*-xy.json"

file_list = glob.glob(pattern)


def write_table(columns, col_names):
    """
    Write an ascii table of columns (sequence of sequences), using col_names as the header
    """
    table = "# " + "\t".join(col_names) + "\n"
    for row in zip(*columns):
        table += "\t".join(row) + "\n"
    return table
        

# Initialize a list for each column in output table
id_ = []
ra = []
dec = []
d = []
pa_star = []
pa_out = []
pa_in = []
r_out = []
r_in = []

for file_name in file_list:
    with open(file_name) as f:
        data = json.load(f)

    # Add this object's data to each output column
    id_.append(data["star"]["id"])
    ra.append(data["star"]["RA"])
    dec.append(data["star"]["Dec"])
    d.append(str(data["star"]["D"]))
    pa_star.append(str(data["star"]["PA"]))
    pa_out.append(str(data["outer"]["PA0"]))
    try:
        pa_in.append(str(data["inner"]["PA0"]))
    except:
        pa_in.append("-"+"\t")
    r_out.append(str(data["outer"]["R0"]))
    try:
        r_in.append(str(data["inner"]["R0"]))
    except:
        r_in.append("-"+"\t")
     
    # D = ???
    # PA_star =
    # PA_in = 
    # PA_out = 
    # R_in = 
    # R_out =

 
# Write output table to a file
with open("arcs-summary.tab", "w") as f:
    f.write(write_table( 
        [id_, ra, dec, d, pa_star, pa_out, pa_in, r_out, r_in], 
        ["Object", "RA", "Dec", "D", "PA_star", "PA_out", "PA_in", "R_out", "R_in"]
    ))
    
