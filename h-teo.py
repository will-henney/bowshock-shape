from __future__ import print_function
from equation6 import Shell
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import bisect
import argparse

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
#4: Plot everything

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

cmd_args=parser.parse_args()
beta = cmd_args.beta
innertype = cmd_args.innertype

#1:
Nth = 800
th_lim = theta_lim(beta)
theta = np.linspace(0,th_lim,Nth)
shell = Shell(beta=beta,innertype=innertype)



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
R_ext,theta=R_ext[R_ext>0], theta[R_ext>0] #avoid bad points (usually the last ones)
#R_circ = R_ext[theta==0]*np.sqrt( 1+( 2*a*( 1-np.cos(theta_c) ) ) / (1-a)**2 )
x_e,y_e = R_ext/R0*np.cos(theta),R_ext/R0*np.sin(theta)
alfa = alpha(x_e,y_e) #slope of CD
alpha_a = -np.arctan(alfa) #Angle of tangent line, defined as positive
psi = theta+alpha_a-0.5*np.pi
M = np.array([2.0,3.0,3.5,4.0]) #Mach Number
theta1 = theta_1(R_ext,theta)


plt.rc("text",usetex=True)
plt.rc("font",family="serif")
f =plt.figure()
ax1 = f.add_subplot(2,1,1)
ax2= f.add_subplot(2,1,2)


for m in M:
    H0 = 3./(4*m**2+1.) #Shell width at axis
    n = 3.4*(4*m**2+1)/(4*m**2+35)
    H = H0*np.cos(psi)**(-n)
    sigma = (np.sin(theta+theta1)*csc(theta)*csc(theta1)*(beta*(1-np.cos(theta))+0.01*m*(1-np.cos(theta1)))**2)/np.sqrt(
         (beta*(theta-np.sin(theta)*np.cos(theta))+(theta1-np.sin(theta1)*np.cos(theta1)))**2+(beta*(np.sin(theta))**2-np.sin(theta1)**2)**2)
    if cmd_args.shell:
        R_in = R_ext/R0 - H/np.cos(psi)
        R_in[R_in<0] = np.nan 
        x_i,y_i = R_in*np.cos(theta),R_in*np.sin(theta)    
        plt.plot(x_i,y_i,"--",label="M={}".format(m))
    else:
        ax1.plot(np.degrees(theta),H,label="M={}".format(m))
        ax2.plot(np.degrees(theta),sigma,label="M={}".format(m))


if cmd_args.shell:
    plt.plot(x_e,y_e,label=r"type={},$\beta$={}".format(innertype,beta))
    plt.grid()
    plt.legend(loc="best",prop=dict(size="x-small"))
    plt.xlabel(r"x/$R_0$")
    plt.xlim(-4.8,1.2)
    plt.ylim(0,6)
    plt.ylabel(r"y/$R_0$")
    plt.title("Shell shape")
    plt.gcf().set_size_inches(8, 8)
    plt.savefig("shell-shape.pdf")

else:
    ax1.legend(loc="best",prop=dict(size="x-small"))
    ax2.legend(loc="best",prop=dict(size="x-small"))
    ax1.grid()
    ax2.grid()
    ax1.set_xlabel(r"$\theta({}^{\circ})$")
    ax2.set_xlabel(r"$\theta({}^{\circ})$")
    ax1.set_ylabel("H")
    ax2.set_ylabel(r"$\sigma/\sigma_0$")
#    ax1.xlim(0,np.degrees(th_lim))
#    ax2.xlim(0,np.degrees(th_lim))
#    ax1.ylim(0,1)
#    ax1.ylim(0,1)
    ax1.set_title(r"H. $\beta={}$,{}".format(beta,innertype))
    ax2.set_title(r"Mass Surface Density".format(beta,innertype))
    f.savefig("H,s-vs-theta.pdf")
