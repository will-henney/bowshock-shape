regionfile = "LV-OIII-positions.reg" 


def extract_data(line):
    coordstring, paramstring = line.split("#")
    shape, numbers = coordstring.split(")")[0].split("(")
    ra, dec = numbers.split(",")[0:2]
    if "tag" in line:
        text = paramstring.split("tag={")[1].split("}")[0]
    elif "text" in line:
        text = paramstring.split("text={")[1].split("}")[0]
    else:
        text = "NO TEXT"
    return ra, dec, text


Shapes = {}

with open(regionfile) as f:
    lines = f.readlines()
    for line in lines: 
        skipthisline = line.startswith("#") \
            or line.startswith("global") \
            or not "#" in line 
        if skipthisline: 
            continue
        ra, dec, text = extract_data(line)
        source = text.split()[0]
        if source not in Shapes:
            Shapes[source] = []
        Shapes[source].append((ra, dec))


for source in Shapes:
    print "******"
    print "***", source, "***"
    print "******"
    print Shapes[source]

