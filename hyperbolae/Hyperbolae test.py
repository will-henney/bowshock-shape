
# In[ ]:

import numpy as np


# In[22]:

pylab.rcParams['figure.figsize'] = (16.0, 8.0)


# In[1]:

def yhyperbol(x, th_inf=45.0):
    tan_th = np.tan(np.radians(th_inf))
    return (np.sqrt(1.0 + x**2*tan_th**2) - 1.0)/tan_th**2

def yparabol(x):
    return 0.5*x**2

def ycircle(x):
    return 1.0 - np.sqrt(1.0 - x**2)


# In[28]:

x = np.linspace(-2.5, 2.5, 5000)


# In[29]:

plot(x, ycircle(x), label="circle")
plot(x, yparabol(x), label="parabola")
for th in 15, 30, 45, 60, 75:
    plot(x, yhyperbol(x, th), label="hyperbola {}".format(th))
axis("equal")
grid(alpha=0.3)
xlabel("x")
ylabel("y")
legend(loc="upper center")


# Out[29]:

#     <matplotlib.legend.Legend at 0x109881fd0>

# image file:

# In[ ]:



