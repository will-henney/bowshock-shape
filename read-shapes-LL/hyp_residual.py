import numpy as np
from equation6 import Shell
import conic_utils
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import argparse
"""
The goal is to find the residual (the relative difference) between a CRW shape and 
the best fit using the points until an angle theta_max (45,60,90,100,120 degrees) and with no
theta_restriction
"""

parser = argparse.ArgumentParser(description="Do the fit with the inner points or the outer points of the bowshock")
parser.add_argument("--out",action="store_true",help="fit with the head of the bowshock or the wings")
cmd_args = parser.parse_args()

beta = [0.002,0.005,0.01,0.05]
theta = np.linspace(0,np.pi,300)
innertype = ["isotropic","proplyd"]
color = ["r","g","b","m","c","k"]
for inn in innertype:
    for b in beta:
        shell = Shell(beta=b,innertype =inn)
        R = shell.radius(theta)
        R,theta = R[R>0],theta[R>0] #filtering bad points
        R0 = R[0]
    #    w = np.zeros_like(R)
    #    w[:-1] = np.diff(np.log(R))/np.diff(np.log(theta)) # may cause problems and it is not neccesary
    #    w[-1]=w[-2]
    #    tan_a= (1+w*np.tan(theta))/(w-np.tan(theta))
    #    thinf = np.arctan(tan_a[-1])
        #completing the bowshock and ordering from -theta_max to theta_max
        R_com,theta_com = np.array([R,R]),np.array([theta,-theta])
        R_com,theta_com = R_com.reshape(2*len(R),),theta_com.reshape(2*len(theta),)
        order = theta_com.argsort()
        R_com,theta_com = R_com[order],theta_com[order]
        #setting the differents upper limits for theta for fitting
        if cmd_args.out:
            theta_mask = [np.radians(45.0),np.radians(60.0),np.radians(90.0),np.radians(100.0),np.radians(120.0),0.0]
        else:
            theta_mask = [np.radians(45.0),np.radians(60.0),np.radians(90.0),np.radians(100.0),np.radians(120.0),np.max(theta_com)]
        for th,col in zip(theta_mask,color):
            if cmd_args.out:
                m = np.abs(theta_com)>= th
            else:
                m = np.abs(theta_com)<= th
            Rm,thm=R_com[m],theta_com[m]
            xfit,yfit = Rm*np.cos(thm),Rm*np.sin(thm)
            #fitting
            Rh,thh,PAh,xh,yh = conic_utils.fit_conic(xfit,yfit,Rh=R0,thh=-30.0,PAh=90.0,xxh=0.0,yyh=0.0)
            xx,yy= conic_utils.world_hyperbola(Rh,thh,PAh,xh,yh)
            #change to R,theta
            RR,tt = np.hypot(xx,yy),np.arctan2(yy,xx)
            #interpolate R from theta
            R_int = interp1d(-tt,RR)#The order in tt is from max_positive to max_negative and theta_com is in the opposite order
            int_mask = np.abs(theta_com)<np.max(tt)
            #calculate the difference between fit and original curve (residual)
            epsilon = np.abs((R_com[int_mask] - R_int(theta_com[int_mask]))/R_com[int_mask])
            #plot residual vs theta for each beta value
            plt.plot(np.degrees(theta_com[int_mask]),epsilon,col,label=str(np.degrees(th)))
        plt.legend(loc="best",prop=dict(size="x-small"))
        plt.grid()
        plt.rc("text",usetex=True)
        plt.rc("font",family="serif")
        plt.xlabel(r"$\theta$ (${}^{\circ}$)")
        plt.ylabel(r"$\epsilon$")
        plt.title(r"Residual for different hyperbola fitting, {} case and $\beta$ = {}".format(inn,b))
        if cmd_args.out:
            plt.savefig("out-{}-beta-{}-res.pdf".format(inn,b))
        else:
            plt.savefig("{}-beta-{}-res.pdf".format(inn,b))
        plt.clf()
