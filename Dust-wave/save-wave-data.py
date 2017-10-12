import sys
import json
import numpy as np
import scanf

import bow_projection as bp
import bow_diagnostic
import dragoid_shape
import ancantoid_shape
import standing_wave


# Parse command line arguments
try:
    BASE_SHAPE_ID = sys.argv[1]
    AMPLITUDE = float(sys.argv[2])
    WAVENUMBER = float(sys.argv[3])
except:
    sys.exit(f"Usage: {sys.argv[0]} BASE_SHAPE_ID AMPLITUDE WAVENUMBER")

# Sensible defaults for the poly fitting to get R_c and R_90
bp.N_NEIGHBORHOOD = 50
bp.DEGREE_POLY_NEIGHBORHOOD = 2
bp.SCALE_NEIGHBORHOOD = 0.03
bp.DEGREE_POLY_NEIGHBORHOOD_90 = 2
bp.SCALE_NEIGHBORHOOD_90 = 0.01


# Choose which base shape according to command-line argument and parse
# out the shape parameters if any
if BASE_SHAPE_ID == "paraboloid":
    base_shape = bp.paraboloid_R_theta
    shape_label = "Paraboloid"
elif BASE_SHAPE_ID == "wilkinoid":
    base_shape = bp.wilkinoid_R_theta
    shape_label = "Wilkinoid"
elif BASE_SHAPE_ID.startswith("cantoid"):
    ibeta, = scanf.scanf("cantoid-beta%d", BASE_SHAPE_ID)
    beta = ibeta / 100000
    base_shape = bp.Spline_R_theta_from_function(
        ngrid=1000, shape_func=bp.cantoid_R_theta, shape_func_pars=(beta,))
    shape_label = rf"Cantoid $\beta = {beta}$"
elif BASE_SHAPE_ID.startswith("ancantoid"):
    ixi, ibeta = scanf.scanf("ancantoid-xi%d-beta%d", BASE_SHAPE_ID)
    xi, beta = ixi / 100, ibeta / 100000
    base_shape = ancantoid_shape.Ancantoid(xi=xi, beta=beta, n=301)
    shape_label = rf"Ancantoid $\xi = {xi:.1f}$, $\beta = {beta}$"
elif BASE_SHAPE_ID.startswith("dragoid"):
    ialpha, = scanf.scanf("dragoid-alpha%d", BASE_SHAPE_ID)
    alpha = ialpha / 100
    base_shape = dragoid_shape.Dragoid(alpha=alpha)
    shape_label = rf"Dragoid $\alpha_\mathrm{{drag}} = {alpha:.2f}$"
    bp.SCALE_NEIGHBORHOOD = 0.2

# Create perturbed shape
shape = standing_wave.StandingWave(base_shape,
                                   amplitude=AMPLITUDE,
                                   wavenumber=WAVENUMBER)


# Create arrays of inclinations and phases
nphase, ninc = 21, 50
phases = np.linspace(0.0, 0.5, nphase)
# Uniform distribution in sin inc
inclinations = np.arcsin(np.linspace(0.0, 1.0, ninc))
# Allocate arrays to hold the (R_c, R_90) results
Rc = np.empty((ninc, nphase))
R90 = np.empty((ninc, nphase))
# Loop over inclinations and phases and save the radii
for i, inclination in enumerate(inclinations):
    for j, phase in enumerate(phases):
        shape.phase = phase
        radii = bp.characteristic_radii_projected(inclination, shape)
        Rc[i, j] = radii['tilde R_c prime']
        R90[i, j] = radii['tilde R_90 prime']


# Save the results
shape_id = f"{BASE_SHAPE_ID}-wave-A{int(100*AMPLITUDE):03d}-N{int(10*WAVENUMBER):02d}"

# Save Rc, R90 arrays in binary format
# To retrieve:
#
# >>> arrays = np.load(FILE.npz)
# >>> Rc = arrays['Rc']
#
# Etcetera
np.savez(shape_id, Rc=Rc, R90=R90)

# Save metadata as JSON 
metadata = {
    'shape_id': shape_id,
    'shape_label': shape_label,
    'wave_label': f"Perturbation $A = {AMPLITUDE}$, $N = {WAVENUMBER}$",
    'ninc': ninc,
    'nphase': nphase,
    'inclinations': list(np.round(np.degrees(inclinations), 1)),
    'phases': list(np.round(phases, 4)),
}
jsonfile = shape_id + ".json"
with open(jsonfile, 'w') as f:
    json.dump(metadata, f, indent=4)

print(jsonfile, end='')
