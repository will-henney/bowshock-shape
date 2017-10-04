import sys
import json
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
    s = np.sign(x0 - 1.0)
    return b*np.sqrt(1.0 + s*((x - x0)/a)**2)

alldata = {}

fit = LevMarLSQFitter()

figfile = sys.argv[0].replace('.py', '.pdf')
sns.set_style('white')
sns.set_color_codes()

fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(6, 4))
efig, eaxes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(6, 4))

alphas = [1.0/4.0, 1.0/2.0, 1.0, 2.0]
results = {'alpha': [0.0] + alphas, 'Rc': [2.0], 'R90': [2.0]}
for alpha, ax, eax in zip(alphas, axes.flat, eaxes.flat):
    astring = f'-alpha{int(100*alpha):03d}.tab'
    alpha_label = fr"$\alpha_\mathrm{{drag}} =  {alpha:.02f}$"
    fitdata = {}
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
    th_head = np.arctan2(y_head, x_head)
    R_head = np.hypot(x_head, y_head)

    # Save the parameters for the head fit
    fitdata['head'] = {'Rc': Rc,
                       'R90': R90,
                       'T': head_conic.b_a}

    # Now find a tail-fit conic
    # We only fit the tail between xswitch and xfar
    xswitch, xfar = -1.0, -8.0
    mtail = (x_bow < xswitch) & (x_bow > xfar)
    mfar = (x_bow <= xfar) & (x_bow > -40.0)

    # Try a more direct approach: fit hyperbola with LM
    model = conic_y_x()
    best_model = fit(model, x_bow[mtail], y_bow[mtail])
    r0_tail = best_model.a.value + best_model.x0.value

    x_tail = np.linspace(r0_tail, -30.0, 1000)
    y_tail = best_model(x_tail)
    th_tail = np.arctan2(y_tail, x_tail)
    R_tail = np.hypot(x_tail, y_tail)

    # Save the parameters for the tail fit
    fitdata['tail'] = {'x0': best_model.x0.value,
                       'a': best_model.a.value,
                       'b': best_model.b.value,
                       'r0': best_model.a.value + best_model.x0.value,
                       'T': best_model.b.value/best_model.a.value}

    # Finally, a third fit to the far tail with a hyperbola
    # model2 = conic_y_x(x0=15.0, a=30.0, b=1.0)
    model2 = conic_y_x(x0=-15.0, a=30.0, b=y_bow[mfar].max())
    best_model2 = fit(model2, x_bow[mfar], y_bow[mfar]) 
    x_far = np.linspace(2.0, -100.0, 1000)
    y_far = best_model2(x_far)
    th_far = np.arctan2(y_far, x_far)
    R_far = np.hypot(x_far, y_far)

    # Save the parameters for the far tail fit
    fitdata['far'] = {'x0': best_model2.x0.value,
                      'a': best_model2.a.value,
                      'b': best_model2.b.value,
                      'r0': best_model2.a.value + best_model2.x0.value,
                      'T': best_model2.b.value/best_model2.a.value}

    # Stash the fit data in the big dict
    alldata[alpha] = fitdata

    # Plot the bow and the two fits
    ax.axvspan(xfar, xswitch, color='k', alpha=0.05)
    ax.plot(x_bow, y_bow, lw=7, alpha=0.3, label='_nolegend_')
    ax.plot(x_far, y_far, ls='-.', color='m', label="Far tail fit")
    ax.plot(x_tail, y_tail, ls=':', lw=2.5, color='r', label="Tail fit")
    ax.plot(x_head, y_head, ls='--', color='orange', label="Head fit")
    ax.text(-14, 0.8, alpha_label, fontsize='small')

    ax.set_aspect('equal', adjustable='box-forced')

    # And plot the errors
    f_R_tail = interp1d(th_tail, R_tail, bounds_error=False)
    e_tail = (f_R_tail(theta) - R_bow)/R_bow
    f_R_head = interp1d(th_head, R_head, bounds_error=False)
    e_head = (f_R_head(theta) - R_bow)/R_bow
    f_R_far = interp1d(th_far, R_far, bounds_error=False)
    e_far = (f_R_far(theta) - R_bow)/R_bow

    # Find angle that corresponds to x = -1
    th1 = interp1d(x_bow, theta)(xswitch)
    # Find angle that corresponds to x = -8
    th2 = interp1d(x_bow, theta)(xfar)
    mh = theta <= th1
    mt = (theta > th1) & (theta <= th2)
    mtt = theta > th2

    eax.axhline(0.0, lw=3, alpha=0.5, color='b')
    eax.axhspan(-0.01, 0.01, color='b', alpha=0.1, ec='none')
    eax.axvspan(np.degrees(th1), np.degrees(th2), color='k', alpha=0.05)

    # Plot each error curve twice, faintly over the bad zone …
    eax.plot(np.degrees(theta[mt | mtt]), e_head[mt | mtt], label='_nolegend_',
             ls='--', color='orange', alpha=0.3)
    # And strongly over the range it should be fitting
    eax.plot(np.degrees(theta[mh]), e_head[mh], label="Head fit",
             ls='--', color='orange', alpha=1.0)

    # Same for the tail fit, but 3 times since we have a bad zone each side
    eax.plot(np.degrees(theta[mh]), e_tail[mh], label='_nolegend_',
             ls=':', lw=2.5, color='r', alpha=0.3)
    eax.plot(np.degrees(theta[mtt]), e_tail[mtt], label='_nolegend_',
             ls=':', lw=2.5, color='r', alpha=0.3)
    eax.plot(np.degrees(theta[mt]), e_tail[mt], label="Tail fit",
             ls=':', lw=2.5, color='r', alpha=1.0)

    # And now for the far-tail fit as well
    eax.plot(np.degrees(theta[mt]), e_far[mt], label='_nolegend_',
             ls='-.', color='m', alpha=0.3)
    eax.plot(np.degrees(theta[mtt]), e_far[mtt], label="Far tail fit",
             ls='-.', color='m', alpha=1.0)

    eax.text(10.0, -0.04, alpha_label, fontsize='small')


for ax in axes[-1, :]:
    ax.set(xlim=[-15, 2.5], xlabel='$x/R_0$')
for ax in axes[:, 0]:
    ax.set(ylim=[0, 8], ylabel='$y/R_0$')
axes[-1,-1].legend(fontsize='small')

for eax in eaxes[-1, :]:
    eax.set(xlim=[0.0, 180.0],
            xlabel=r"Polar angle: $\theta$, degrees",
            xticks=[0, 30, 60, 90, 120, 150, 180])
for eax in eaxes[:, 0]:
    eax.set(ylim=[-0.05, 0.05],
            ylabel=r"Fractional error: $\delta R / R$")
eaxes[-1,-1].legend(fontsize='small')
eaxes[-1, 1].text(10.0, -0.9/100, r"$|\delta R/R| < 1\%$", 
                  color='b', fontsize='x-small')

# print(results)
#     ax.plot(theta, R, label=fr"${alpha:.2f}$")

# ax.legend(title=r"$\alpha_\mathrm{drag}$")
# ax.axhline(2.0, color='k', alpha=0.3, lw=1)
# ax.axvline(90.0, color='k', alpha=0.3, lw=1)
# ax.set(xlim=[0.0, 180.0], ylim=[0.0, 14.0],
#        xlabel=r'$\theta$', ylabel=r'$R$')

sns.despine(fig)
fig.tight_layout()
fig.savefig(figfile)

sns.despine(efig)
efig.tight_layout()
efig.savefig(figfile.replace('.pdf', '-error.pdf'))

jsonfile = figfile.replace('-Rc-R90.pdf', '-fitdata.json')
with open(jsonfile, 'w') as f:
    json.dump(alldata, f, indent=4)

print(figfile, end='')
