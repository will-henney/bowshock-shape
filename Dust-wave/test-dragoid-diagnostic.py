import sys
import numpy as np
import bow_projection as bp
import dragoid_shape
import bow_diagnostic

try:
    alpha = float(sys.argv[1])
except:
    sys.error(f"Usage: {sys.argv[0]} ALPHA")

shape = dragoid_shape.Dragoid(alpha=alpha)

th_inf = bp.theta_infinity(shape)
inclinations = np.linspace(0.0, th_inf - np.pi/2, 30)
Rc, R90 = bow_diagnostic.Rcp_R90p(inclinations, shape)

result = [['inc', 'R_c', 'R_90'], None]
result += list(zip(np.degrees(inclinations).astype(int), Rc, R90))
print(result)
