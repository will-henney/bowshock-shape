import json

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

shelldata = json.load(open("rc-r0.json"))

for modeldata in shelldata.items():
    beta,param = modeldata
    r0 = np.array(param["R0'"])
    rc = np.array(param["Rc"])
    inc = np.array(param["inc"])
    plt.plot(np.array(float(beta)), r0[inc==0.0]/rc[inc==0.0], "r.")

def y(x, b, d): 
    c = np.exp(-1.0/b)
    return y0*(np.exp(-((x**d)/b)) - c)/(1.0 - c)

y0 = 0.66
b = 50.0
d = 0.49
x = np.linspace(0.0, 0.1, 500)
plt.plot(x, y(x, b, d))
plt.grid()
plt.xlabel("beta")
plt.ylabel("1/A")
# plt.xlim(0.0, 1.0)
# plt.ylim(0.0, 1.0)
plt.savefig("AVSb.pdf")
print "b = ", b
print "d = ", d
print "exp(-1/b) = ", np.exp(-1.0/b)
print "y0^-1 = ", 1./y0
