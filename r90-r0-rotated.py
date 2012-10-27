import numpy as np
import matplotlib.pyplot as plt
from equation6 import Shell
from scipy.optimize import fsolve,bisect
"""
In this script, we'll create the proyection of proplyd in the sky plane, assuming the plane of proplyd is rotated respect plane of sky by an angle i
"""
def theta_lim(beta): 
    "Asymptotic opening angle from CRW eq (28)"
    def f(theta):
	"Function to be zeroed: f(theta) = 0 for theta = theta_lim"
	return theta - np.tan(theta) - np.pi/(1.0 - beta)
    return bisect(f, 0.5*np.pi+0.01, np.pi)

def omega(r,t):
    """
    Will's version of omega(theta), looks like work better.
    uses d(log(r))/dtheta = 1/r * dr/dtheta (ingenious)
    """
    womega = np.zeros_like(r)
    womega[:-1] = np.diff(np.log(r))/np.diff(theta)
    return womega

def R90_finder(x,y):
    """
    This fuction finds a plausible value for R90, without calculing Theta_perp. In the ratated frame, r90 is y such that x = 0. So i find this value
    looking for the place where x has a sign change, if this function does not find any sign change, returns a NaN. If x and y don't have the same size, 
    returns error (for completness) 
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
            out = 0.5*(y[i]+y[i+1])
            return out
    return np.nan



#First, design a set of beta values to calculate R(theta) (Well, first than that, desgign a theta array)
#parser = argparse.ArgumentParser(description="Rotation angle")#
#parser.add_argument("--iangle", "-i", default=0.0, type=float, help="Rotation of z axis in radians")
#cmdargs = parser.parse_args()
#i = cmdargs.iangle

beta = [0.3, 0.1, 0.01, 0.001, 1.e-4, 1.e-5, 1.e-7]
Nth = 200
Ninc = 50
#beta = np.linspace(0.01,0.5,20)

innertype = ['isotropic','proplyd']
# innertype = ['isotropic']


ls = dict(isotropic = '-', proplyd = '--')
for inn in innertype:
    print "******{} case******".format(inn)
    for b in beta:
        print "beta = ", b
        thlim = theta_lim(b)
        theta = np.linspace(0,thlim,Nth)
        i = np.linspace(-0.98*(thlim - 0.5*np.pi),0.98*(thlim - 0.5*np.pi),Ninc)
        R0 = list()
        R90 = list()
        shell = Shell(beta=b, innertype=inn)    
        R = shell.radius(theta)
        w = omega(R,theta)
        for j in i:
    #       x,y,z = R*np.cos(THETA),R*np.sin(THETA)*np.cos(PHI),R*np.sin(THETA)*np.sin(PHI)
            SenodePhiT = np.tan(j) * ( ( 1+ w*np.tan(theta) )/( w-np.tan(theta) ) )
            SenodePhiT[np.abs(SenodePhiT)>=1.] =np.nan #other way to set mask
    #       mask = (np.sin(PHI-0.5*DPHI) <= SenodePhiT)  &  (SenodePhiT <= np.sin(PHI +0.5*DPHI))
    #       xi,yi,zi = x*np.cos(j*np.pi/180.)-z*np.sin(j*np.pi/180.),y,x*np.sin(j*np.pi/180.)+z*np.cos(j*np.pi/180.)Convention sign in article Henney et al
            xi = R*np.cos(theta)*np.cos(j)-R*np.sin(theta)*SenodePhiT*np.sin(j)
            yi = R*np.sin(theta)*np.sqrt(1-SenodePhiT**2)
            mask = np.isfinite(xi) & np.isfinite(yi)
            xim, yim=xi[mask], yi[mask] #Removing nan elements from xi and yi
            R0.append(xim[0]) # Looks like is a secure criterium ( At least in the beta range (0,0.1] ). For high i and high beta (~0.5), odd things happen
            R90.append(R90_finder(xim,yim))
        if inn == "isotropic":
            label = 'beta={}'.format(b)
        else:
            label = None
        plt.plot(R0,R90, ls[inn], label = label)


plt.legend(loc="upper left")
#    plt.axis([0,0.5,0,1.5])
plt.xlabel("R_0 / D")
plt.ylim(0.0, 4.0)
plt.ylabel("R_90 / D")
plt.title("Perpendicular versus parallel bowshock radii")
plt.savefig("combined-shell-test2-R0-R90.pdf")
plt.clf()
