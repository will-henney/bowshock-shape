import numpy as np
import argparse
regionfile = "LV-OIII-positions.reg" 


def extract_data(line):
    coordstring, paramstring = line.split("#")
    shape, numbers = coordstring.split(")")[0].split("(")
    ra, dec = numbers.split(",")[0:2]
    if "tag" in line:
        text = paramstring.split("tag={")[1].split("}")[0]
        label = paramstring.split("tag={")[1].split()[1].split("}")[0]
    elif "text" in line:
        text = paramstring.split("text={")[1].split("}")[0]
        label = paramstring.split("text={")[1].split('-')[-1].split("}")[0]
    else:
        text = "NO TEXT"
        label = 'None'
    return ra, dec, text,label

def extract_centers(shape):

    """
    Extracting the position of a desired proplyd
    Arguments: Shape: Array of tuples with (ra,dec,label) as each element
    """
    for tup in shape:
        if tup[-1] == ' best':
            ra,dec = tup[0],tup[1]
            break
        elif tup[-1] == ' best position':
            ra,dec = tup[0],tup[1]
            break
        else:
            ra,dec = np.nan,np.nan
    return ra,dec,tup[-1]
    
def extract_shocks(shape):
    """
    Extracting the coordinates of the whole shock for a specific proplyd
    Arguments: Same for the extracts_centers function
    """
    r,d = [],[]
    for tup in shape:
        if tup[-1] == 'shock':
            r.append(tup[0])
            d.append(tup[1])
        else:
            continue
    return np.array(r),np.array(d)

parser = argparse.ArgumentParser(description="Name of desired proplyd(s) to measure radii")
parser.add_argument('--prop','-p',default = 'LV4',type = str,help = 'Name of proplyd or proplyds to measure radii')
cmdargs = parser.parse_args()

Shapes = {}

with open(regionfile) as f:
    lines = f.readlines()
    for line in lines: 
        skipthisline = line.startswith("#") \
            or line.startswith("global") \
            or not "#" in line 
        if skipthisline: 
            continue
        ra, dec, text,label = extract_data(line)
        shr,smin,ssec = ra.split(':')
        hr,mn,sec = float(shr),float(smin),float(ssec)
        ra_arcsec = 15*(3600.*hr+60.*mn+sec)
        sdeg,samin,sasec = dec.split(':')
        deg,amin,asec = float(sdeg),float(samin),float(sasec)
        dec_arcsec = 3600.*deg + np.sign(deg)*(60*amin + asec)
        source = text.split()[0]
        if source not in Shapes:
            Shapes[source] = []
        Shapes[source].append((ra_arcsec, dec_arcsec,label))

#Extract individually the coords from th1C

th1c_ra,th1c_dec = Shapes['th1C'][0][0],Shapes['th1C'][0][1]


for LVs in cmdargs.prop.split('-'):
    ra_cen,dec_cen,lab = extract_centers(Shapes[LVs])
    ra_sck,dec_sck = extract_shocks(Shapes[LVs])                           # Extract coordinates from proplyd and shock
    d= np.sqrt( (ra_cen - th1c_ra)**2 + (dec_cen-th1c_dec)**2 )            # distance between proplyd and star
    tp = np.arctan( np.abs( (dec_cen-th1c_dec)/(ra_cen-th1c_ra) ) )        #Angle between LOS star-proplyd and x axis
    tss=   np.arctan( np.abs( (dec_sck-th1c_dec)/(ra_sck-th1c_ra) ) ) - tp #angle between LOS star-shock and LOS star-proplyd
    rt = np.sqrt( (dec_cen-dec_sck)**2 + (ra_cen-ra_sck)**2 )              #distance between star and shock
    rs = np.sqrt( ( ra_sck - th1c_ra )**2 + (dec_sck - th1c_dec )**2 )     #distance between shock and star
    sin_ts = np.sin(tss)*rs/rt                                             #Sines Law
    ts = np.degrees( np.arcsin(sin_ts)  )
    print LVs,'rt/d=',rt.max()/d,'theta=',ts[rt==rt.max()]
    print LVs,'r0/d=',rt.min()/d,'theta=',ts[rt==rt.min()]


