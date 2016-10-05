import numpy as np
from astropy.modeling.models import custom_model
from astropy.modeling.fitting import LevMarLSQFitter

@custom_model
def hyperbola_ratio(theta, a=1.0, x0=1.0, tau=1.0, D=2.0):
    """Ratio th1/th for a hyperbola with center `x0`, scale `a`, and
    opening angle `tau` = tan theta_infty.  Complementary angle th1 is
    measured from point D along the x axis

    """
    # First we find the parametric variable t
    tau_cot_theta = tau/np.tan(theta)
    cosht = ((x0 / a)
             * (1 - tau_cot_theta*np.sqrt(1 + (a**2/x0**2)*(tau_cot_theta**2 - 1)))
             / (1 - tau_cot_theta**2))
    # Now we can find x, y
    x = x0 - a*cosht
    y = a*tau*np.sqrt(cosht**2 - 1.0)
    # Then the th1 angle
    theta1 = np.arctan2(y, D - x)
    # And finally return the ratio
    return theta1/theta
