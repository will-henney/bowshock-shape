"Utility functions for CRW96 bowshocks"
import numpy as np
import scipy
from scipy.optimize import fsolve, bisect
from numpy import sin, cos, tan, sqrt, abs, pi

def cot(t):
    "Cotangent"
    return 1./tan(t)

def csc(t):
    "Cosecant"
    return 1./sin(t)


def t1(t,b): 
    "Approximate solution for theta_1 in CRW"
    return sqrt(7.5*(-1.0+sqrt(1.0+0.8*b*(1.0-t/tan(t)))))

def R(t,b):
    "Radius of CRW bowshock"
    if t == 0.0:
	return 1.0
    elif b == 0.0: 
	return Ri(t)
    else:
	return sin(t1(t,b))/sin(t + t1(t,b))/R0(b)

def R0(b):
    "On-axis distance or CRW bowshock"
    return sqrt(b)/(1.0+sqrt(b))

def Ri(t):
    "Radius of Wilkin bowshock"
    return sqrt(3*(1.0-t/tan(t)))/sin(t)

def theta_lim(beta): 
    "Asymptotic opening angle from CRW eq (28)"
    def f(theta):
	"Function to be zeroed: f(theta) = 0 for theta = theta_lim"
	return theta - tan(theta) - pi/(1.0 - beta)
    return bisect(f, 0.5*pi+0.01, pi)
 
def omega(t,b):
    "R^{-1} d R/d theta"
#     numerator = b*(cot(t)-t*csc(t)**2)
#     denominator = 1 + sin(t1(t,b))*cot(t)
#     denominator *= cos(t1(t,b)) - t1(t,b)*csc(t1(t,b))
#     return numerator/denominator - cot(t + t1(t,b))
    return cot(t1(t,b))*dt1dt(t,b) -  cot(t + t1(t,b))*(1 + dt1dt(t,b))

def vshell(t,b,alpha):
    """Velocity along the bowshock shell normalized to that of the inner wind, v_w.
    Extra parameter alpha is ratio of wind velocities: alpha = v_w/v_w1 """
    st = sin(t) ; ct = cos(t); st1 = sin(t1(t,b)); ct1 = cos(t1(t,b))
    top = sqrt( 
	(b*(t - st*ct) + (t1(t,b) - st1*ct1))**2 + (b*st**2 - st1**2)**2
	)
    bottom = 2*(b*(1-ct) + alpha*(1-ct1))
    return top/bottom

def dt1dt(t,b):
    def f(t):
	return cot(t) - t/sin(t)**2
    return b*f(t)/f(t1(t,b))


# 
# Functions for the projected axis line: y' = 0
# 
# def xax(t, b, inc):
#     "x'(y' = 0) returns 2-tuple since there are two solutions"

# def vax_los(t,b,alpha, inc):
#     """
#     Line-of-sight velocity of gas at positions along projected axis (y' = 0)
#     Returns a tuple since the LOS crosses the shell in two points
#     """
#     vs = vshell(t,b,alpha)
    
#     return ( , )


# 
# Functions for the tangent line
# 
def sphit(t,b,inc):
    "sin phi_t"
    return tan(inc)*(1 + omega(t,b)*tan(t))/(omega(t,b) - tan(t))

def xt(t, b, inc):
    "x'_t"
    return R(t,b)*(cos(t)*cos(inc)-sin(t)*sphit(t,b,inc)*sin(inc))

def yt(t, b, inc):
    "y'_t"
    return R(t,b)*sin(t)*sqrt(1-sphit(t,b,inc)**2)

def vt_los(t,b,alpha, inc):
    "Line-of-sight velocity of gas at positions along tangent line"
    om = omega(t, b)
    return vshell(t,b,alpha) * sqrt(1 + om**2) * sin(inc) / (sin(t) - om*cos(t)) 

#
# Particular points on the tangent line:
#                                        par is at y'_t = 0
#                                        perp is at x'_t = 0
#
def theta_par_approx(b, inc):
    "Small angle approximation to thpar"
    return inc / (1.0 - 0.4*(1.0 + 1.5*sqrt(b)))

def theta_par(b, inc, thmin = 1.e-4, thmax=None, thlim=0.5*pi):
    """Minimum theta, corresponding to y'_t = 0
    """
    if (inc == 0.0): return 0.0
    if thmax == None:
	# calculate thmax if it is not passed as argument
	thmax = theta_lim(b)
    tani = abs(tan(inc))
    def f(theta):
	"Function to be zeroed: f(theta) = 0 for theta = theta_par"
	om = omega(theta,b)
	return sin(theta)*(1.0 - om*tani) - cos(theta)*(om + tani)
    print '%.2e %.4f f(%.2e) is %.2e' % (b, inc, thmin, f(thmin))
    print '%.2e %.4f f(%.2e) is %.2e' % (b, inc, thmax, f(thmax))
    while f(thmin)*f(thmax) > 0.0:
	if 1.1*thmax < thlim:
	    thmax = 1.1*thmax	# hunt upwards
	else:
	    thmax = 0.99*thlim
	    thmin = 0.9*thmin	# hunt downwards
	print '%.2e %.4f f(%.2e) is %.2e' % (b, inc, thmax, f(thmax))
    print thmin, thmax, f(thmin), f(thmax)
    return bisect(f, thmin, thmax)

def theta_perp(b, inc, thmin=None, thmax=None):
    "theta corresponding to x'_t = 0"
    if (inc == 0.0): return 0.5*pi
    if thmax == None:
	# calculate this is not passed as argument
	thmax = theta_lim(b)
    if thmin == None:
	thmin = theta_par(b, inc, thmax)
    sinsq2i = sin(2*inc)**2
    cossqi = cos(inc)**2
    def f(theta):
	"Function to be zeroed: f(theta) = 0 for theta = theta_perp"
	om = omega(theta,b)
	return cot(theta) - (1 - sqrt(1 + om**2 * sinsq2i))/(2*om*cossqi)
    if f(thmin)*f(thmax) > 0.0:
	thmin = 0.5*pi
    print "theta_perp: ", thmin, thmax, f(thmin), f(thmax)
    return bisect(f, thmin, thmax)

def Rpar(b, inc, thpar=None):
    "Apparent radius along projected symmetry axis"
    if thpar == None:
	thpar = theta_par(b, inc)
    return R(thpar, b)*cos(inc - thpar)

def Rperp(b, inc, thperp=None):
    "Apparent radius perpendicular to projected symmetry axis"
    if thperp == None:
	thperp = theta_perp(b, inc)
    return R(thperp, b)*sin(thperp)*sqrt(1.0 - sphit(thperp,b,inc)**2)

