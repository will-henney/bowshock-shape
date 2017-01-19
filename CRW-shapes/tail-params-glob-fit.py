import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from astropy.table import Table

# ************* Read the astropy table with the data ************************

tab = Table.read("conic-head-tail-fit-finegrid.tab",format="ascii.tab")


# ************* Fit the residual for x0 and fit a - x0 **********************

def x0_trend(beta):
    return 0.7*beta**-0.55

# Function to fit y= x0/x0_trend: y = a*log(beta)**2 + b*log(beta) + c
# Function to fit z= a-x0: z = d + k*log(beta)

# *************  for beta < 5e-4 didn't get a good fit ***********************

m1 = tab["beta"] > 5e-5

b = tab["beta"][m1]
y = tab["x0"][m1]/x0_trend(b)
z = tab["x0"][m1]-tab["a"][m1]
logb = np.log10(b)

# ---------> Till here, everything looks fine

xis = set(tab["xi"])
nxi = len(xis) 
# ***** Separating isotropic case from the anisotropic case ********

xi = np.array(sorted(xis))

iso_case = xi[0]

xi_aniso = xi[1:]

# ********** Fit the anisotropic parameters x0/x0_trend and x0 - a ***************************
#*********** and Plot the analytic x0 -a and x0 against the data *****************************
C1_x0, C2_x0,C3_x0,C4_x0 = [],[],[],[]
C1_a,C2_a,C3_a = [],[],[]

# *************** Set figures styles **************************
fig, (ax1,ax2) = plt.subplots(2,1,sharex=True)
FIG, (AX1,AX2) = plt.subplots(2,1,sharex=True)
sns.set_style("whitegrid")
colors = sns.color_palette("Blues_d",n_colors=nxi)

for x,c in zip(xi_aniso,colors):
    Label = r"$\xi={:.2f}$".format(float(x))
    m2 = tab["xi"][m1] == x
    p_x0 = np.polyfit(logb[m2],y[m2],3)
    p_xa = np.polyfit(logb[m2],z[m2],2)
    c1_x0,c2_x0,c3_x0,c4_x0 = p_x0[0],p_x0[1],p_x0[2],p_x0[-1]
    c1_a,c2_a,c3_a = p_xa[0],p_xa[1],p_xa[-1]
    C1_x0.append(c1_x0)
    C2_x0.append(c2_x0)
    C3_x0.append(c3_x0)
    C4_x0.append(c4_x0)
    C1_a.append(c1_a)
    C2_a.append(c2_a)
    C3_a.append(c3_a)
    x0_analytic = c1_x0*logb[m2]**3 + c2_x0*logb[m2]**2 + c3_x0*logb[m2]+c4_x0
    x0a_analytic = c1_a*logb[m2]**2 + c2_a*logb[m2]+c3_a
    AX1.plot(b[m2],y[m2],c=c,alpha=0.5)
    AX1.plot(b[m2],x0_analytic,c=c,linestyle="--",label=Label)
    AX2.plot(b[m2],z[m2],alpha=0.5,c=c)
    AX2.plot(b[m2],x0a_analytic,c=c,linestyle="--",label=Label)
    
#AX1.legend()
#AX2.legend()
AX1.set(xscale="log",yscale="linear",ylabel=r"$x_0$")
AX2.set(xscale="log",yscale="linear",xlabel=r"$\beta$",ylabel=r"$x_0-a$")
FIG.savefig("fit-tail-params-1.pdf")

# ***** The isotropic case requires a sigle fit ******


m_iso = tab["xi"] == iso_case
p0_iso = np.polyfit(logb[m_iso],y[m_iso],3)
pa_iso = np.polyfit(logb[m_iso],z[m_iso],2)

xi_an_to_fl = xi_aniso.astype(np.float)




# ****************  Find the dependence of the fit parameters with xi ************************

# Note: All the parameters looks like they have linear dependencies with xi
#F0_C1 = np.polyfit(xi_an_to_fl,C1_x0,1)
#F0_C2 = np.polyfit(xi_an_to_fl,C2_x0,1)
#F0_C3 = np.polyfit(xi_an_to_fl,C3_x0,1)

#Fa_C1 = np.polyfit(xi_an_to_fl,C1_a,1)
#Fa_C2 = np.polyfit(xi_an_to_fl,C2_a,1)

#F0_C1_fit = F0_C1[0]*xi_an_to_fl + F0_C1[-1]
#F0_C2_fit = F0_C2[0]*xi_an_to_fl + F0_C2[-1]
#F0_C3_fit = F0_C3[0]*xi_an_to_fl + F0_C3[-1]

#Fa_C1_fit = Fa_C1[0]*xi_an_to_fl + Fa_C1[-1]
#Fa_C2_fit = Fa_C2[0]*xi_an_to_fl + Fa_C2[-1]
# ****** Plot the parameters and the respective fits********

#ax1.plot(xi_an_to_fl,C1_x0,c="r",label=None,alpha=0.5)
#ax1.plot(xi_an_to_fl,C2_x0,c="g",label=None,alpha=0.5)
#ax1.plot(xi_an_to_fl,C3_x0,c="b",label=None,alpha=0.5)

#ax1.plot(1.0,p0_iso[0],"ro",label=r"Isotropic Fit: $C_1 = {:.4f}$".format(p0_iso[0]))
#ax1.plot(1.0,p0_iso[1],"go",label=r"Isotropic Fit: $C_2 = {:.4f}$".format(p0_iso[1]))
#ax1.plot(1.0,p0_iso[-1],"bo",label=r"Isotropic Fit: $C_3 = {:.4f}$".format(p0_iso[-1]))

#ax1.plot(xi_an_to_fl,F0_C1_fit,"r--",
#         label=r"Fit: $C_1 = {:.4f}\xi + {:.4f}$".format(F0_C1[0],F0_C1[-1]))
#ax1.plot(xi_an_to_fl,F0_C2_fit,"g--",
#         label=r"Fit: $C_2 = {:.4f}\xi + {:.4f}$".format(F0_C2[0],F0_C2[-1]))
#ax1.plot(xi_an_to_fl,F0_C3_fit,"b--",
#         label=r"Fit: $C_3 = {:.4f}\xi + {:.4f}$".format(F0_C3[0],F0_C3[-1]))

#ax1.set(xscale="linear",xlim=[0.2,1.05],yscale="linear",ylim=[None,None])
#ax1.legend(loc="best",title=r"$\frac{x_0}{0.7\beta^{-0.55}}=C_1\log_{10}\beta^2+C_2\log_{10}\beta+C_3$",
#           fontsize="xx-small",ncol=2)

#ax2.plot(xi_an_to_fl,C1_a,c="b",label=None,alpha=0.5)
#ax2.plot(xi_an_to_fl,C2_a,c="r",label=None,alpha=0.5)
#ax2.plot(1.0,pa_iso[0],"bo",label=r"Isotropic Fit: $C_1 = {:.4f}$".format(pa_iso[0]))
#ax2.plot(1.0,pa_iso[-1],"ro",label=r"Isotropic Fit: $C_2 = {:.4f}$".format(pa_iso[-1]))

#ax2.plot(xi_an_to_fl,Fa_C1_fit,"b--",
#         label=r"Fit: $C_1 = {:.4f}\xi + {:.4f}$".format(Fa_C1[0],Fa_C1[-1]))
#ax2.plot(xi_an_to_fl,Fa_C2_fit,"r--",
#         label=r"Fit: $C_2 = {:.4f}\xi + {:.4f}$".format(Fa_C2[0],Fa_C2[-1]))


#ax2.set(xscale="linear",xlim=[0.2,1.1],yscale="linear",ylim=[None,None],xlabel=r"$\xi$")
#ax2.legend(loc="best",title=r"$a-x_0=C_1\log\beta+C_2$",fontsize="x-small")

#fig.set_size_inches(4,6)
#fig.tight_layout()
#fig.savefig("tail_params_fit.pdf")

# -------> Till here, everything works, except for a few warnings I don't understand





