import sys
sys.path.append("projected")
import bowfuncs as bow
import numpy as np
import matplotlib.pyplot as plt


N = 100
beta = np.linspace(0.01,0.5,20)
i = np.linspace(0,20,9)
for j in i:
    R0 = list()
    R90 = list()
    for b in beta:
        R0.append(bow.Rpar(b,j))
        R90.append(bow.Rperp(b,j))
    plt.plot(R0,R90,'o-',label = 'i={}'.format(j))
plt.legend(loc="upper left")
#plt.axis([0,0.5,0,1.5])
plt.xlabel("R_0 / D")
plt.ylabel("R_90 / D")
plt.title("Perpendicular versus parallel bowshock radii")
plt.savefig("will-shell-test-R0-R90.png")
plt.clf()    
