from __future__ import print_function
from equation6 import Shell
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import bisect
import argparse
import json

"""
Main goal: Trace the theoretical inner and outer shell in bowshock
Hypothesis:
1.-  v*cos(\psi)/c \equiv M_\perp.
 Where v is the inner wind velocity, 
c is the sound speed, beta_a is the angle between the velocity vector and the normal
 vector in the shell and M_\perp is the mach number perpendicular to the shell.
2.- h \propto (\cos\psi)**{-n}.
Where h is the shell width and n varies with Mach Number
So, we conclude that h/h_0 = (cos(\psi))**{-n} where h_0 is the shell width at theta=0
3.- R_in = R - h/\cos\psi
"""

#Steps:
#1: Create a theta array and create a bowshock shell using CRW formalism
#2: Calculate \psi for each theta: we can use geometrical arguments
#3: Calculate h and R_in
#4: Change to the plane of sky reference frame and get h'
#5: Plot everything

def csc(t):
    """
    Gives the cosecant of an angle
    """
    return 1./np.sin(t)

def theta_1(R,th):
    """
    Theta_1 is the angle between the axis and a point in the 
    bowshock measured from the massive star, it's needed to calculate the 
    superficial density from CRW formalism
    """
    return np.arctan2(R*np.sin(th),1.-R*np.cos(th))

def theta_lim(beta):
    "Asymptotic opening angle from CRW eq (28)"
    def f(theta):
        "Function to be zeroed: f(theta) = 0 for theta = theta_lim"
        return theta - np.tan(theta) - np.pi/(1.0 - beta)
    return bisect(f, 0.5*np.pi+0.01, np.pi)

def omega(r, t):
    """
    Will's version of omega(theta), looks like work better.  uses
    d(log(r))/dtheta = 1/r * dr/dtheta (ingenious)
    """
    womega = np.zeros_like(r)
    womega[:-1] = np.diff(np.log(r))/np.diff(t)
    return womega

def alpha(x,y):
    """
    tangent of alpha is the drivative of y respect x
    """
    a = np.zeros_like(x)
    a[:-1] = np.diff(y)/np.diff(x)
    a[-1] = a[-2]
    return a



parser = argparse.ArgumentParser(description= """Inputs: wind type and inner wind parameter "beta" """)
parser.add_argument("--beta",type=float,default = 0.01,help= "winds momentum-rate ratio")
parser.add_argument("--innertype",type=str,choices=("isotropic","proplyd"),default="proplyd",help="inner wind model")
parser.add_argument("--shell",action="store_true",help="If true, plot the shell shape, else plot h vs theta")
parser.add_argument("--inc",type=float,default = 15.0,help="inlination angle (in degrees) between proplyd reference frame and plane of sky reference frame")
cmd_args=parser.parse_args()
beta = cmd_args.beta
innertype = cmd_args.innertype

#1:
Nth = 800
th_lim = theta_lim(beta)
theta = np.linspace(0,np.radians(130),Nth)
shell = Shell(beta=beta,innertype=innertype)
inc = np.radians(cmd_args.inc)


#Add the analytic fit found for A = Rc/R0
y0 = {"proplyd":0.66,"isotropic": 0.585}
b = {"proplyd":50.0,"isotropic":5.0}
c = np.exp(-1./b[innertype])
d = {"proplyd": 0.49, "isotropic": 0.5}
A =  (1-c)/(y0[innertype]*(np.exp(-beta**d[innertype]/b[innertype])-c))
a = (A-1)/A


#The respective radii (CRW and circular)
R_ext = shell.radius(theta)
R0=R_ext[0]
#R_ext,theta=R_ext[R_ext>0], theta[R_ext>0] #avoid bad points (usually the last ones)
RN=R_ext/R0
#R_circ = R_ext[theta==0]*np.sqrt( 1+( 2*a*( 1-np.cos(theta_c) ) ) / (1-a)**2 )
x_e,y_e = RN*np.cos(theta),RN*np.sin(theta)
alfa = alpha(x_e,y_e) #slope of CD
alpha_a = -np.arctan(alfa) #Angle of tangent line, defined as positive

phitan= alfa*np.tan(inc) #Tangent line is at an azimutal angle \phi_t, this is \sin\phi_t
x_pe,y_pe = (RN)*(np.cos(theta)*np.cos(inc) - np.sin(theta)*np.sin(inc)*phitan),(RN)*(np.sin(theta)*np.sqrt(1-phitan**2))
#x_pe,y_pe = x_pe[np.isfinite(x_pe) & np.isfinite(y_pe)],y_pe[np.isfinite(x_pe) & np.isfinite(y_pe)]
R_pe = np.hypot(x_pe,y_pe)
theta_p = np.arctan2(y_pe,x_pe)
alpha_pa = -np.arctan(alpha(x_pe,y_pe))
psi_p = theta_p+alpha_pa-0.5*np.pi

psi = theta+alpha_a-0.5*np.pi
M = np.array([2.0,3.0,3.5,4.0]) #Mach Number
cm = ["r","g","b","m"]
theta1 = theta_1(R_ext,theta)


plt.rc("text",usetex=True)
plt.rc("font",family="serif")

out = {"innertype":innertype,
"beta":beta,
"theta":list(theta),
"inc":cmd_args.inc,
"psi":list(psi),
"CD Radius":list(RN),
"indexes":[3.4*(4*m**2+1)/(4*m**2+19) for m in M],
"H0":[3./(4*m**2+1.) for m in M],
"Rc/R0":A
}

for i,m in enumerate(M):
    H0 = 3./(4*m**2+7.) #Shell width at axis
    n = 3.4*(4*m**2+1)/(4*m**2+19)
    H = H0*np.cos(psi)**(-n)
    R_in = RN - H/np.cos(psi)
    R_in[R_in<=0] = np.nan 
    x_i,y_i = R_in*np.cos(theta),R_in*np.sin(theta)
    alpha_in = alpha(x_i,y_i)
    phitan_i = np.tan(inc)*alpha_in
    x_pi,y_pi = R_in*(np.cos(theta)*np.cos(inc) - np.sin(theta)*np.sin(inc)*phitan_i),R_in*(np.sin(theta)*np.sqrt(1-phitan_i**2))
    #mask=np.isfinite(x_pi) & np.isfinite(y_pi)
    #x_pi,y_pi = x_pi[mask],y_pe[mask]
    R_pi = np.hypot(x_pi,y_pi)
    H_p = (R_pe-R_pi)*np.cos(psi_p)
    out["H'(M={})".format(m)] = list(H_p)
    if cmd_args.shell:    
        plt.plot(x_i,y_i,cm[i]+"--",label="M={}".format(m))
	plt.plot(x_pi,y_pi,cm[i]+":",alpha=0.5,label="M={} projected, i={}".format(m,cmd_args.inc))
    else:
        plt.plot(np.degrees(theta),H,cm[i],label="M={}".format(m))
        plt.plot(np.degrees(theta_p),H_p,cm[i]+"--",label="M={} projected, i={}".format(m,cmd_args.inc))


if cmd_args.shell:
    plt.plot(x_e,y_e,"k-",lw=2,label=r"type={}, $\beta$={}".format(innertype,beta))
    plt.plot(x_pe,y_pe,"k-",lw=2,alpha =0.5,label="projected, i={}".format(cmd_args.inc))
    plt.grid()
    plt.legend(loc="best",prop=dict(size="x-small"))
    plt.xlabel(r"x'/$R_0$")
    plt.xlim(-4,1.2)
    plt.ylim(0,4)
    plt.ylabel(r"y'/$R_0$")
    plt.title("Shell shape")
    plt.gcf().set_size_inches(8, 8)
    plt.savefig("shell-shape-{}-{}-inc-{}.pdf".format(innertype,beta,cmd_args.inc))

else:
    plt.grid()
    plt.legend(loc="best",prop=dict(size="x-small"))
    plt.xlabel(r"$\theta'({}^{\circ})$")
    plt.ylabel("H'")
    plt.xlim(0,np.degrees(th_lim))
    plt.ylim(0,1.0)
    plt.title(r"H'. $\beta={}$, {}".format(beta,innertype))
    plt.savefig("H-{}-{}-inc-{}.pdf".format(beta,innertype,cmd_args.inc))

#save the output in a json file



jsonfile = "H-{}-{}-inc-{}.json".format(beta,innertype,cmd_args.inc)

with open(jsonfile,"w") as f:
    json.dump(out,f,indent=2)
