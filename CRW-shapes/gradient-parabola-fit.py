
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

# ################## t parameter for parabola ####################
#Q2 = 1
#Q1 = 2/np.tan(theta)
#Q0 = -2*x0p/Rcp
#tp = (-Q1 + np.sqrt(Q1**2 - 4*Q2*Q0))/(2*Q2)
##################################################################

################ t parameter for ellipse ########################
#C2 = a**2 + b**2/np.tan(theta)**2
#C1 = 2*b*(a-1)/np.tan(theta)
#C0 = 1-2*a
#sin_te = (-C1 + np.sqrt(C1**2 - 4*C2*C0))/(2*C2)
#cos_te = np.sqrt(1. - sin_te**2)
################################################################# 

############### Plot gradient and linear fit ####################
x = R*np.cos(theta)
y = R*np.sin(theta)
dxdy = np.diff(x)/np.diff(y)
y_min = R_point(np.radians(90))*np.sin(np.radians(90)) 
y_max = R_point(np.radians(150))*np.sin(np.radians(90))
mfit, y_ref = par_fit(dxdy,y[:-1],y_min,y_max)
dxdy_line = mfit*y[:-1] + y_ref
f = plt.figure()
#ax1 = f.add_subplot(3, 1, 1, adjustable="box", aspect=1)
ax1 = f.add_subplot(1, 1, 1, adjustable="box", aspect=0.1)
ax1.plot(y[:-1],dxdy,label="Wilkin")
ax1.plot(y[:-1],dxdy_line,label="Linear fit")
ax1.legend()
ax1.set_xlabel(r"$y$")
ax1.set_ylabel(r"$\frac{dx}{dy}$")
f.savefig("gradient-par-test.pdf")
#################################################################
