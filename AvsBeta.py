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
    plt.plot(np.array(float(beta)), rc[inc==0.0]/r0[inc==0.0], "r.")
plt.grid()
plt.xlabel("beta")
plt.ylabel("A")
plt.savefig("AVSb.pdf")
