
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
def wilkin_tail_params(t1,t2):
    def pointR(t):
        return np.sqrt(3*(1-t/np.tan(t))/np.sin(t)**2)
    x1 = pointR(t1)*np.cos(t1)
    x2 = pointR(t2)*np.cos(t2)
    y1 = pointR(t1)*np.sin(t1)
    y2 = pointR(t2)*np.sin(t2)
    Rc = -0.5*(y1**2-y2**2)/(x1-x2)
    x0 = 0.5*(x1+x2) + 0.25*(y1**2+y2**2)/Rc
    return Rc,x0

#def R(T):
#    if t == 0:
theta1 = np.radians(110)
theta2 = np.radians(150)
Rcp,x0p = wilkin_tail_params(theta1,theta2)
theta = np.linspace(0,np.pi,500,endpoint=False)
R = np.sqrt(3*(1-theta/np.tan(theta))/np.sin(theta)**2)
x = R*np.cos(theta)
y = R*np.sin(theta)

# ################## t parameter for parabola ####################
Q2 = 1
Q1 = 2/np.tan(theta)
Q0 = -2*x0p/Rcp
tp = (-Q1 + np.sqrt(Q1**2 - 4*Q2*Q0))/(2*Q2)
##################################################################

xp = -0.5*Rcp*tp**2 + x0p
yp = Rcp*tp
R90 = np.sqrt(3.)
Rc = 5./3
Tc = 2*Rc - R90**2
a = Rc/Tc
b = Rc/np.sqrt(Tc)

################ t parameter for ellipse ########################
C2 = a**2 + b**2/np.tan(theta)**2
C1 = 2*b*(a-1)/np.tan(theta)
C0 = 1-2*a
sin_te = (-C1 + np.sqrt(C1**2 - 4*C2*C0))/(2*C2)
cos_te = np.sqrt(1. - sin_te**2)
################################################################# 

xe = a*cos_te + 1-a
ye = b*sin_te
f = plt.figure()
ax1 = f.add_subplot(2,1,1,adjustable="box",aspect=1)
ax1.plot(x,y,"k-",lw=2,alpha=0.5)
ax1.plot(xp,yp)
ax1.plot(xe,ye)
ax1.set_xlim(-50,1)
#plt.gca().set_aspect("equal","box")
plt.ylim(0,20)
# Plot residuals

R_par = np.sqrt(xp**2 + yp**2)
Re = np.sqrt(xe**2 + ye**2)
epsilon = np.abs(R - R_par)/R
epsilone = np.abs(R - Re)/R
ax2 = f.add_subplot(2,1,2,adjustable="box",aspect=100)
#ax2.plot(np.degrees(theta),np.degrees(theta_par))
#ax2.plot(np.degrees(theta),R_par)
#ax2.plot(np.degrees(theta),R2,"k-")
#ax2.plot(np.degrees(theta),Re)
ax2.plot(np.degrees(theta),epsilon)
ax2.plot(np.degrees(theta),epsilone)
#ax2.set_xlim(0,160)
ax2.set_ylim(-0.1,1)
f.set_size_inches(6,6)
f.savefig("2-points-par-test.pdf")
