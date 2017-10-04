import sys
import numpy as np
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
import seaborn as sns
from astropy.table import Table
import astropy.modeling.fitting
from astropy.modeling.models import custom_model
from astropy.modeling.fitting import LevMarLSQFitter
sys.path.append('../conic-projection')
from conproj_utils import Conic

@custom_model
def conic_y_x(x, x0=-3.0, a=5.0, b=3.0):
    return b*np.sqrt(1.0 - ((x - x0)/a)**2)

fit = LevMarLSQFitter()

figfile = sys.argv[0].replace('.py', '.pdf')
sns.set_style('white')
sns.set_color_codes()

fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(6, 4))

alphas = [1.0/4.0, 1.0/2.0, 1.0, 2.0]
results = {'alpha': [0.0] + alphas, 'Rc': [2.0], 'R90': [2.0]}
for alpha, ax in zip(alphas, axes.flat):
    astring = f'-alpha{int(100*alpha):03d}.tab'
    t = Table.read('dust-couple-stream' + astring, format='ascii.tab')
    dth = np.pi/len(t)
    theta = t['theta'] + 0.5*dth
    # Mask to select only the near-apex region
    m = np.degrees(theta) <= 30.0
    # Reflect through origin to ensure an even function
    thth = np.concatenate([-theta[m][::-1], theta[m]])
    RR = np.concatenate([t['R'][m][::-1], t['R'][m]])
    # Fit polynomial to find R0 and Rc
    a, b, c = np.polyfit(thth, RR, deg=2)
    assert(b == 0.0, 'Linear term should be non-zero, but is not')
    R0 = c
    Rc = 1.0/(1.0 - 2.0*a/c)

    # Find R90 by linear interpolation
    f = interp1d(theta, t['R'], kind='linear')
    R90 = f(0.5*np.pi)/R0
    results['Rc'].append(Rc)
    results['R90'].append(R90)

    R_bow = t['R']/R0
    x_bow = R_bow*np.cos(theta)
    y_bow = R_bow*np.sin(theta)

    # Define the head-fit conic section that corresponds to (Rc, R90)
    arg = 2*Rc - R90**2
    thc = np.sign(arg)*np.arctan(np.sqrt(np.abs(arg)))
    head_conic = Conic(A=Rc, th_conic=np.degrees(thc))
    t = head_conic.make_t_array()
    x_head = head_conic.x(t)
    y_head = head_conic.y(t)


    # Now find a tail-fit conic
    mtail = (x_bow < -1.0) & (x_bow > -8.0)

    # Try a more direct approach: fit hyperbola with LM
    model = conic_y_x()
    best_model = fit(model, x_bow[mtail], y_bow[mtail])
    print(best_model.a, best_model.b, best_model.x0)
    r0_tail = best_model.a.value + best_model.x0.value

    # # Add the other branch of the tail
    # xx = np.concatenate([x_bow[mtail][::-1], x_bow[mtail]])
    # yy = np.concatenate([-y_bow[mtail][::-1], y_bow[mtail]])
    # # This is reverse polynomial: x(y)
    # tail_coeffs = np.polyfit(yy, xx, deg=2)
    # assert(tail_coeffs[1] == 0.0, tail_coeffs)
    # ptail = np.poly1d(tail_coeffs)

    x_tail = np.linspace(r0_tail, -10.0, 200)
    y_tail = best_model(x_tail)

    ax.plot(x_bow, y_bow, lw=3, alpha=0.5)
    ax.plot(x_head, y_head, ls='--')
    ax.plot(x_tail, y_tail, ls=':', lw=2.5, color='r')
    ax.set_aspect('equal', adjustable='box-forced')

for ax in axes[-1, :]:
    ax.set(xlim=[-8, 2], xlabel='$x/R_0$')

for ax in axes[:, 0]:
    ax.set(ylim=[0, 6], ylabel='$y/R_0$')

# print(results)
#     ax.plot(theta, R, label=fr"${alpha:.2f}$")

# ax.legend(title=r"$\alpha_\mathrm{drag}$")
# ax.axhline(2.0, color='k', alpha=0.3, lw=1)
# ax.axvline(90.0, color='k', alpha=0.3, lw=1)
# ax.set(xlim=[0.0, 180.0], ylim=[0.0, 14.0],
#        xlabel=r'$\theta$', ylabel=r'$R$')

sns.despine()
fig.tight_layout()
fig.savefig(figfile, dpi=300)
print(figfile, end='')
