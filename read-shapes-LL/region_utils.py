region_hdr_lines = [
    "# Region file format: DS9 version 4.1",
    "global color=red dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1",
    "fk5",
]
  
def region_box_to_string(ra="5:34:40.800", dec="-5:22:43.00",
                         width=50, height=50, angle=0.0, 
                         text="LL2", color="orange"):
    return "box({},{},{}\",{}\",{}) # color={} text={{{}}}".format(
        ra, dec, width, height, angle, color, text)

def region_circle_to_string(ra="5:34:40.800", dec="-5:22:43.00",
                            radius=10.0, text="LL2", color="yellow"):
    return "circle({},{},{}\") # color={} text={{{}}}".format(
        ra, dec, radius, color, text)

def region_point_to_string(ra="5:34:40.800", dec="-5:22:43.00",
                            ptype="diamond", text="LL2", color="yellow"):
    return "point({},{}) # point={} color={} text={{{}}}".format(
        ra, dec, ptype, color, text)
# point(5:35:11.391,-5:22:42.90) # point=diamond color=white text={inner}

