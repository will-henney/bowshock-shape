
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
theta1 = np.radians(150)
theta2 = np.radians(179)
Rcp,x0p = wilkin_tail_params(theta1,theta2)
theta = np.linspace(0,np.pi,300,endpoint=False)
R = np.sqrt(3*(1-theta/np.tan(theta))/np.sin(theta)**2)
x = R*np.cos(theta)
y = R*np.sin(theta)
t = np.linspace(0,200,1000)
xp = -0.5*Rcp*t**2 + x0p
yp = Rcp*t
R90 = np.sqrt(3.)
Rc = 5./3
Tc = 2*Rc - R90**2
a = Rc/Tc
b = Rc/np.sqrt(Tc)
xe = a*np.cos(theta) + 1-a
ye = b*np.sin(theta)
f = plt.figure()
ax1 = f.add_subplot(2,1,1)
ax1.plot(x,y,"k-",lw=2,alpha=0.5)
ax1.plot(xp,yp)
ax1.plot(xe,ye)
ax1.set_xlim(-25,x0p)
#plt.gca().set_aspect("equal","box")
#plt.ylim(0,45)
# Plot residuals

R_par = np.sqrt(xp**2+yp**2)
theta_par = np.arctan2(yp,xp)
#Draw again Wilkin bowshock and ellipse using theta_par
R2 = np.sqrt(3*(1-theta_par/np.tan(theta_par))/np.sin(theta_par)**2)
xe2, ye2 = a*np.cos(theta_par) + 1 -a, b*np.sin(theta_par)
Re = np.sqrt(xe2**2+ye2**2)
epsilon = np.abs(R2 - R_par)
epsilone = np.abs(R2 - Re)
ax2 = f.add_subplot(2,1,2)
#ax2.plot(np.degrees(theta),np.degrees(theta_par))
ax2.plot(np.degrees(theta_par),R_par)
ax2.plot(np.degrees(theta_par),R2,"k-")
ax2.plot(np.degrees(theta_par),Re)
#ax2.plot(np.degrees(theta_par),epsilon)
#ax2.plot(np.degrees(theta_par),epsilone)
#ax2.set_xlim(0,160)
ax2.set_ylim(0,100)
f.set_size_inches(6,6)
f.savefig("2-points-par-test.pdf")
