
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def R_point(t):
    """
    Evaluate R at a singular value for theta
    """
    return np.sqrt(3*(1-t/np.tan(t))/np.sin(t)**2)

def par_fit(D,y,y1,y2):
    """
    Do linear fit for dx/dy in the [y1,y2] range
    """
    mask = (y1 < y) & (y < y2) 
    fit = np.polyfit(y[mask],D[mask],1)
    m = fit[0]
    yref = fit[1]
    return m, yref 
theta = np.linspace(0,np.pi,400,endpoint=False)
R = np.sqrt(3*(1-theta/np.tan(theta))/np.sin(theta)**2)


############### Plot gradient and linear fit ####################
x = R*np.cos(theta)
y = R*np.sin(theta)
dxdy = np.diff(x)/np.diff(y)
y_min = R_point(np.radians(90))*np.sin(np.radians(90)) 
y_max = R_point(np.radians(150))*np.sin(np.radians(150))
mfit, y_ref = par_fit(dxdy,y[:-1],y_min,y_max)
dxdy_line = mfit*y[:-1] + y_ref
f = plt.figure()
#ax1 = f.add_subplot(3, 1, 1, adjustable="box", aspect=1)
ax1 = f.add_subplot(3, 1, 1, adjustable="box", aspect=0.1)
ax1.plot(y[:-1], dxdy, label="Wilkin")
ax1.plot(y[:-1], dxdy_line, label="Linear fit")
ax1.legend()
ax1.set_xlabel(r"$y$")
ax1.set_ylabel(r"$\frac{dx}{dy}$")
#################################################################

################## Plot parabola and ellipse ####################
Rcp = -1./mfit
y1 = y_max
x1 = y1/np.tan(np.radians(150))
x0p = x1 + 0.5*y1**2/Rcp
# ################## t parameter for parabola ####################
Q2 = 1
Q1 = 2/np.tan(theta)
Q0 = -2*x0p/Rcp
tp = (-Q1 + np.sqrt(Q1**2 - 4*Q2*Q0))/(2*Q2)
xp = -0.5*Rcp*tp**2 + x0p
yp = Rcp*tp
##################################################################

################ t parameter for ellipse ########################
R90 = np.sqrt(3.)
Rce = 5./3
Tc = 2*Rce - R90**2
a = Rce/Tc
b = Rce/np.sqrt(Tc)
C2 = a**2 + b**2/np.tan(theta)**2
C1 = 2*b*(a-1)/np.tan(theta)
C0 = 1-2*a
sin_te = (-C1 + np.sqrt(C1**2 - 4*C2*C0))/(2*C2)
cos_te = np.sqrt(1. - sin_te**2)
xe = a*cos_te + (1-a)
ye = b*sin_te
################################################################# 
ax2 = f.add_subplot(3, 1, 2, adjustable="box", aspect=1)
ax2.plot(x, y, "k-", lw=2, alpha=0.7, label="Wilkin")
ax2.plot(xp, yp, label="Parabola fit")
ax2.plot(xe, ye, label="Elliptic head")
ax2.legend()
ax2.set_xlabel(r"$x$")
ax2.set_ylabel(r"$y$")  
ax2.set_xlim(-100,1)
ax2.set_ylim(0,20)
#################################################################

##################### Plot residuals ############################
R_par = np.sqrt(xp**2 + yp**2)
Re = np.sqrt(xe**2 + ye**2)
epsilon = np.abs(R - R_par)/R
epsilone = np.abs(R - Re)/R
ax3 = f.add_subplot(3, 1, 3, adjustable="box", aspect=50)
ax3.plot(np.degrees(theta), epsilon, label="Parabolic Tail")
ax3.plot(np.degrees(theta), epsilone, label="Elliptic Head")
ax3.set_ylim(-0.1,1) 
ax3.set_xlabel(r"$\theta$ (deg)")
ax3.set_ylabel(r"$\epsilon$")
#################################################################
f.savefig("gradient-par-test.pdf")
