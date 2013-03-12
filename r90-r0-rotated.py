"""
In this script, we create the projection of proplyd in the sky
plane, assuming the proplyd axis is at an angle i with respect to
plane of sky
"""

import argparse
from StringIO import StringIO

import numpy as np
from scipy.optimize import fsolve,bisect,leastsq
import matplotlib
import matplotlib.pyplot as plt

from equation6 import Shell


#
# Command line arguments 
#
parser = argparse.ArgumentParser(
    description="""Find parallel and perpendicular projected radii for
    bowshock models"""
    )#
parser.add_argument("--observations", type=str, 
                    choices=("jorge_mir", "will_mir","jorge_60","jorge_45","jorge_curv"), default="will_mir", 
                    help="Which version of the observational measurements to use")
parser.add_argument("--yaxis", type=str, 
                    choices=("R90/D", "R90/R0"), default="R90/D", 
                    help="Which quantity to plot on the y axis")
parser.add_argument('--tobs',type = float,default = 90,help='projected bowshock radius at theta = tobs respect to axis')
cmd_args = parser.parse_args()
isNormalized = cmd_args.yaxis == "R90/R0"

def theta_lim(beta): 
    "Asymptotic opening angle from CRW eq (28)"
    def f(theta):
	"Function to be zeroed: f(theta) = 0 for theta = theta_lim"
	return theta - np.tan(theta) - np.pi/(1.0 - beta)
    return bisect(f, 0.5*np.pi+0.01, np.pi)

def omega(r,t):
    """
    Will's version of omega(theta), looks like work better.  uses
    d(log(r))/dtheta = 1/r * dr/dtheta (ingenious)
    """
    womega = np.zeros_like(r)
    womega[:-1] = np.diff(np.log(r))/np.diff(theta)
    return womega

def Rt_finder(x,y,t):
    """
    This function finds a plausible value for R90, without calculating
    Theta_perp. In the rotated frame, r90 is y such that x = 0. So i
    find this value looking for the place where x has a sign change,
    if this function does not find any sign change, returns a NaN. If
    x and y don't have the same size, returns error (for completeness)
    """
    size = np.size(x)
    sizey = np.size(y)
    assert size == sizey, "x and y array sizes must be equal" 

    # # If inclination is < 90, then x decreases monotonically
    # npos = len(x[x > 0.0])         # number of positive x values
    # if npos > 0 and npos > size:
    #     out = 0.5*(y[npos-1]+y[npos])
    #     return out

    for i in range(0,size-1):
        s = (np.degrees(np.arctan2(y[i],x[i])) - t)*( np.degrees( np.arctan2(y[i+1],x[i+1]) ) - t )
        if s<0:
            slope = (y[i+1]-y[i])/(x[i+1]-x[i])
            xout = (y[i]-slope*x[i])/( np.tan( np.radians(t) ) -slope )
            yout = y[i] + slope*(xout-x[i])
            return np.sqrt(xout**2+yout**2)
    return np.nan

def calc_R(xc, yc):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((xim_fit-xc)**2 + (yim_fit-yc)**2)

#@countcalls
def f_2(c):
    """ calculate the algebraic distance between the 2D points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(*c)
    return Ri - Ri.mean()

def wf_2(c):
    """ calculate the algebraic distance between the 2D weighted points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(*c)
    return u*(Ri - Ri.mean())


#First, design a set of beta values to calculate R(theta) (Well, first
#than that, desgign a theta array)




beta = [0.16, 0.08, 0.04, 0.02, 0.01, 0.005, 0.002]
Nth = 400
Ninc = 200

innertype = ['isotropic','proplyd']
# innertype = ['isotropic']

params = {
    "font.family": "serif",
    "text.usetex": True,
    "text.latex.preamble": [r"\usepackage[varg]{txfonts}"],
    "figure.figsize": (5, 5),
    }
matplotlib.rcParams.update(params)



# Jorge's measured radii for proplyds
# LV1		05 35 16.83	-05 23 26.23	1"	2"	6.44"	0.16	0.31
# LV2		05 35 16.74	-05 23 16.51	1.87"	2.6"**	7.68"	0.24	0.3
# LV3		05 35 16.28	-05 23 16.69	2"	3.14"	6.82"	0.29	0.46
# LV4		05 35 16.06	-05 23 24.42	1.08"	1.78"	6.18"	0.17	0.29
# LV5		05 35 15.83	-05 23 22.59	1.94"	2" **	9.46"	0.21	0.21

# Will's measured radii for proplyds
"""
| Object |    D' |  R'_0 | R'_90A | R'_90B | R0/D | R90/D | R90/R0 | Notes       |
|--------+-------+-------+--------+--------+------+-------+--------+-------------|
| LV1    | 6.651 | 0.830 |  2.260 |  2.260 | 0.12 |  0.34 |   2.83 | S side only |
| LV2    | 7.777 | 2.078 |  2.442 |  2.442 | 0.27 |  0.31 |   1.15 | E side only |
| LV3    | 6.977 | 2.057 |  3.422 |  3.422 | 0.29 |  0.49 |   1.69 | E side only |
| LV4    |  6.22 | 1.098 |  1.763 |   2.05 | 0.18 |  0.31 |   1.72 |             |
| LV5    | 9.426 | 1.819 |  2.831 |  3.189 | 0.19 |  0.32 |   1.68 |             |
| LV2b   | 7.024 | 0.815 |  2.146 |  2.146 | 0.12 |  0.31 |   2.58 | W side only |
      #+TBLFM: $6=$3/$2 ;f2::$7=0.5 ($4 + $5)/$2; f2::$8=$7/$6 ; f2
"""

obs_data = dict(
    jorge_mir = """
    1 0.16 0.31 -20 20 
    2 0.24 0.3  20 0 
    3 0.29 0.46 20 20 
    4 0.17 0.29 0 -20 
    5 0.21 0.21 20 -20
    """,
    will_mir = """
    1 0.12 0.34 -20 20 
    2 0.27 0.31  20 0 
    3 0.29 0.49 20 20 
    4 0.18 0.31 -20 -20 
    5 0.19 0.32 20 -20
    20 0.12 0.31 -20 0 
    """,
    jorge_60 = """
    3 0.33 0.43 20 20
    4 0.19 0.22 -20 -20
    5 0.21 0.28 20 -20
    20 0.09 0.09 -20 20
    141 0.06 0.10 20 20
    180 0.05 0.06 20 20
    """,
    jorge_45 = """ 
    2 0.24 0.32 20 20
    3 0.33 0.39 20 20
    4 0.19 0.22 0 20
    5 0.21 0.24 0 -20
    20 0.09 0.11 -20 20
    141 0.06 0.08 0 20
    176 0.13 0.14 -20 0
    """,
    jorge_curv = """
    2 0.33 0.55 20 20
    3 0.24 0.23 20 -20
    4 0.19 0.37 0 -20
    5 0.21 0.33 20 20
    20 0.09 0.23 0 20
    141 0.06 0.12 20 20
    176 0.13 0.16 -20 0
    180 0.05 0.11 -20 20
    """
    )


obs_labels, obs_R0, obs_R90, obs_dx, obs_dy = \
    np.loadtxt(StringIO(obs_data[cmd_args.observations]), unpack=True)

obs_label = ['LV2','LV3','LV4','LV5','LV2b','141-301','176-341','180-331']
lw = dict(isotropic = 2, proplyd = 3)
opacity = dict(isotropic = 0.3, proplyd = 0.7)

colors = "bgrmkcy"
for inn in innertype:
    print "******{} case******".format(inn)
    for b, col in zip(beta, colors[:len(beta)]):
        print "beta = ", b
        thlim = theta_lim(b) # I have to calculate theta_lim for
                             # proplyd case
        theta = np.linspace(0,thlim,Nth)
        if inn == 'proplyd':
            # In proplyd case, try all inclinations
            inclinations = np.linspace(0.0, 0.5*np.pi, Ninc)   
        else:
            inclinations = np.linspace(0.,0.98*(thlim - 0.5*np.pi),Ninc)
            print "Maximum Inclination: ", np.degrees(inclinations[-1])

        # Mask to select inclinations close to multiples of 15 degrees
        every15 = np.abs((np.degrees(inclinations) + 7.5) % 15.0 - 7.5) <= 0.5*np.degrees(inclinations[1])

        # Initialize the data arrays to NaNs
        R0 = np.zeros_like(inclinations)*np.nan
        Rt = np.zeros_like(inclinations)*np.nan
        Rc = np.zeros_like(inclinations)*np.nan

        shell = Shell(beta=b, innertype=inn)    
        R = shell.radius(theta)
        w = omega(R,theta)
        for i, inc in enumerate(inclinations):
            SenodePhiT = np.tan(inc) * ( ( 1+ w*np.tan(theta) )/( w-np.tan(theta) ) )
            SenodePhiT[np.abs(SenodePhiT)>=1.] =np.nan #other way to set mask

            # Correct for fact that observed radii are normalised by
            # D' = D cos(inc)
            xi = (R/np.cos(inc))*(np.cos(theta)*np.cos(inc) 
                                - np.sin(theta)*SenodePhiT*np.sin(inc))
            yi = (R/np.cos(inc))*np.sin(theta)*np.sqrt(1-SenodePhiT**2)
            mask = np.isfinite(xi) & np.isfinite(yi)
            xim, yim=xi[mask], yi[mask] #Removing nan elements from xi
                                        #and yi
            xim_fit,yim_fit = xim[xim>=0],yim[xim>=0]
            c0 = 0,0 # initail guess for the circle fit to calculate curvature radius
            try:
                c_fit,ier = leastsq(f_2,c0)
                xfit,yfit = c_fit
                Ri_2 = calc_R(xfit, yfit)
            except TypeError:
                print "Maximum Inclination:",np.degrees(inc)
                break
                #skip invalid inputs for fitting and skip to next beta value
            try:
                # Looks like is a secure criterion ( At least in the
                # beta range (0,0.1] ). For high i and high beta
                # (~0.5), odd things happen
                R0[i] = xim[0] 
                Rc[i] = Ri_2.mean() # The fit uses ALL the bowshock data
            except IndexError:
                print "Maximum Inclination: ", np.degrees(inc)
                break        # ignore inclinations with no valid
                             # solution and skip remaining incs
 
        label = r'\(\beta={}\)'.format(b) if inn == "proplyd" else None
        Y = Rc/R0 if isNormalized else Rc
            
        # First, plot a line with all the inclinations
        plt.plot(R0, Y, '-', linewidth=lw[inn], c=col, label=label, alpha=opacity[inn])
        # Second, add symbols every 15 degrees
        print np.degrees(inclinations[every15])
        plt.plot(R0[every15], Y[every15], '.', c=col, alpha=opacity[inn])
        

# Add the observations to the plot
obs_Y = obs_R90/obs_R0 if isNormalized else obs_R90
plt.plot(obs_R0, obs_Y, "ko")
for label, x, y, dx, dy in zip(obs_label, obs_R0, obs_Y, obs_dx, obs_dy):
    ha = "right" if dx < 0 else "left"
    va = "top" if dy < 0 else "bottom"
    plt.annotate(
        "{}".format(label), 
        xy = (x, y), xytext = (dx, dy),
        textcoords = 'offset points', ha = ha, va = va,
        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0')
        )


plt.xlabel(r"\(R'_0 / D'\)")
# avoid numbering the origin
epsilon = 1.e-6
if isNormalized:
    plt.xlim(0.0 + epsilon, 0.5)
    plt.ylim(0.0 - epsilon, 8 + epsilon)
    plt.ylabel(r"\(R'_{c} / R'_0\)")
    plt.legend(loc="best", ncol=2, prop=dict(size="x-small"))
else:
    plt.xlim(0.0 + epsilon, 0.35)
    plt.ylim(0.0 + epsilon, 1.0)
    plt.ylabel(r"\(R'_{c} / D'\)")
    plt.legend(loc="best")

suffix = "-norm" if isNormalized else ""

plt.title("Curvature radii versus parallel bowshock radii")
plt.savefig("combined-shell-R0-Rc{}.pdf".format(suffix))
plt.clf()
