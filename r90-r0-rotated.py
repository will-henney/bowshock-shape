import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from equation6 import Shell
from StringIO import StringIO
from scipy.optimize import fsolve,bisect
"""
In this script, we create the projection of proplyd in the sky
plane, assuming the proplyd axis is at an angle i with respect to
plane of sky
"""
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

def R90_finder(x,y):
    """
    This fuction finds a plausible value for R90, without calculating
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
        s = x[i]*x[i+1]
        if s<0:
            out = y[i] - x[i]*( ( y[i+1]-y[i] )/ (x[i+1]-x[i]) ) #Linear interpolation
            return out
    return np.nan




#First, design a set of beta values to calculate R(theta) (Well, first
#than that, desgign a theta array)
#parser = argparse.ArgumentParser(description="Rotation angle")#
#parser.add_argument("--iangle", "-i", default=0.0, type=float, help="Rotation of z axis in radians")
#cmdargs = parser.parse_args()
#i = cmdargs.iangle



beta = [0.16, 0.08, 0.04, 0.02, 0.01, 0.005]
Nth = 400
Ninc = 200

innertype = ['isotropic','proplyd']
# innertype = ['isotropic']

params = {
    "font.family": "serif",
    "text.usetex": True,
    "text.latex.preamble": [r"\usepackage[varg]{txfonts}"],
    "figure.figsize": (5, 10),
    }
matplotlib.rcParams.update(params)



# Measured radii for proplyds
# LV1		05 35 16.83	-05 23 26.23	1"	2"	6.44"	0.16	0.31
# LV2		05 35 16.74	-05 23 16.51	1.87"	2.6"**	7.68"	0.24	0.3
# LV3		05 35 16.28	-05 23 16.69	2"	3.14"	6.82"	0.29	0.46
# LV4		05 35 16.06	-05 23 24.42	1.08"	1.78"	6.18"	0.17	0.29
# LV5		05 35 15.83	-05 23 22.59	1.94"	2" **	9.46"	0.21	0.21

obs_data = StringIO("""
    1 0.16 0.31 -20 20 
    2 0.24 0.3  20 0 
    3 0.29 0.46 20 20 
    4 0.17 0.29 0 -20 
    5 0.21 0.21 20 -20
    """)
obs_labels, obs_R0, obs_R90, obs_dx, obs_dy = np.loadtxt(obs_data, unpack=True)

lw = dict(isotropic = 2, proplyd = 3)
opacity = dict(isotropic = 0.3, proplyd = 0.7)

colors = "bgrmky"
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
        R90 = np.zeros_like(inclinations)*np.nan

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
            try:
                # Looks like is a secure criterion ( At least in the
                # beta range (0,0.1] ). For high i and high beta
                # (~0.5), odd things happen
                R0[i] = xim[0] 
                R90[i] = R90_finder(xim,yim)
            except IndexError:
                print "Maximum Inclination: ", np.degrees(inc)
                break        # ignore inclinations with no valid
                             # solution and skip remaining incs
 
        label = r'\(\beta={}\)'.format(b) if inn == "proplyd" else None
        # First, plot a line with all the inclinations
        plt.plot(R0,R90, '-', linewidth=lw[inn], c=col, label=label, alpha=opacity[inn])
        # Second, add symbols every 15 degrees
        print np.degrees(inclinations[every15])
        plt.plot(R0[every15], R90[every15], '.', c=col, alpha=opacity[inn])
        

# Add the observations to the plot
plt.plot(obs_R0, obs_R90, "ko")
for label, x, y, dx, dy in zip(obs_labels, obs_R0, obs_R90, obs_dx, obs_dy):
    ha = "right" if dx < 0 else "left"
    va = "top" if dy < 0 else "bottom"
    plt.annotate(
        "LV{:0.0f}".format(label), 
        xy = (x, y), xytext = (dx, dy),
        textcoords = 'offset points', ha = ha, va = va,
        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))


plt.legend(loc="upper left")
#    plt.axis([0,0.5,0,1.5])
plt.xlabel(r"\(R'_0 / D'\)")
# avoid numbering the origin
epsilon = 1.e-6
plt.xlim(0.0 + epsilon, 0.35)
plt.ylim(0.0 + epsilon, 1.0)
plt.ylabel(r"\(R'_{90} / D'\)")
plt.title("Perpendicular versus parallel bowshock radii")
plt.savefig("combined-shell-test2-R0-R90.pdf")
plt.clf()
