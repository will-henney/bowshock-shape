from __future__ import print_function
import json

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

shelldata = json.load(open("rc-r0.json"))
shelltype = ["isotropic","proplyd"]
shellcolor= ["r.","g."]

for inn,col in zip(shelltype,shellcolor):
    label = inn
    for modeldata in shelldata[inn].items():
        beta,param = modeldata
        r0 = np.array(param["R0'"])
        rc = np.array(param["Rc"])
        inc = np.array(param["inc"])
        plt.plot(np.array(float(beta)), rc[inc==0.0]/r0[inc==0.0], col,label=label)
        label =None

def y(x, b, d): 
    c = np.exp(-1.0/b)
    return y0*(np.exp(-((x**d)/b)) - c)/(1.0 - c)
print("Proplyd parameters")
y0 = 0.66
b = 50.0
d = 0.49
x = np.linspace(0.0, 0.1, 500)
plt.plot(x, 1.0/y(x, b, d),label = "proplyd fit")
print ("b = ", b                      )
print ("d = ", d                      )
print ("exp(-1/b) = ", np.exp(-1.0/b) )
print ("y0^-1 = ", 1./y0              )

print ("Isotropic parameters")
y0 = 0.585
b = 5
d = 0.5
plt.plot(x, 1.0/y(x, b, d),"g-",label="isotropic fit")

plt.plot(x, 1.5/(1.0 - np.sqrt(x)), "r-", label='analytic')
plt.legend(loc="best")
plt.grid()
plt.xlabel("beta")
plt.ylabel("A")
# plt.xlim(0.0, 1.0)
# plt.ylim(0.0, 1.0)
plt.savefig("AVSb.pdf")
print ("b = ", b                      )
print ("d = ", d                      )
print ("exp(-1/b) = ", np.exp(-1.0/b) )
print ("y0^-1 = ", 1./y0              )  
