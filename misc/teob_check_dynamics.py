import numpy as np
import json
from nrtools.utils.utils import nu_to_q, q_to_nu, csv_reader
from watpy.utils import num
import EOBRun_module
import matplotlib
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import random
matplotlib.rcParams['text.usetex']= True
matplotlib.rcParams['font.serif']= 'Palatino' 
matplotlib.rcParams['font.size']= 15 #28

from matplotlib.collections import LineCollection

# Phase of the signal
def phase(signal):
    return -np.unwrap(np.angle(signal))

# omega of the signal
def phi_dot(time, signal):
    phase_to_diff = phase(signal)
    return num.diff1(time, phase_to_diff, pad=True)

def multiline(xs, ys, c, ax=None, **kwargs):
    """Plot lines with different colorings

    Parameters
    ----------
    xs : iterable container of x coordinates
    ys : iterable container of y coordinates
    c : iterable container of numbers mapped to colormap
    ax (optional): Axes to plot on.
    kwargs (optional): passed to LineCollection

    Notes:
        len(xs) == len(ys) == len(c) is the number of line segments
        len(xs[i]) == len(ys[i]) is the number of points for each line (indexed by i)

    Returns
    -------
    lc : LineCollection instance.
    """

    # find axes
    ax = plt.gca() if ax is None else ax

    # create LineCollection
    segments = [np.column_stack([x, y]) for x, y in zip(xs, ys)]
    lc = LineCollection(segments, **kwargs)

    # set coloring of line segments
    #    Note: I get an error if I pass c as a list here... not sure why.
    lc.set_array(np.asarray(c))

    # add lines to axes and rescale 
    #    Note: adding a collection doesn't autoscalee xlim/ylim
    ax.add_collection(lc)
    ax.autoscale()
    return lc

def modes_to_k(modes):
    """
    Map multipolar (l,m) -> linear index k
    """
    return [int(x[0]*(x[0]-1)/2 + x[1]-2) for x in modes]


modes = [[2,1],[2,2],[3,2],[3,3],[4,4]]
k = modes_to_k(modes)

mode = str(modes_to_k([[2,2]])[0])

q = 2
mns = 1.4
mbh = q*mns
M = mbh+mns
e = 0.2

pars = {
                'M'                  : M,
                'q'                  : q,
                'chi1x'              : 0,
                'chi1y'              : 0,
                'chi1z'              : 0.7,
                'chi2'               : 0.,
                'LambdaAl2'          : 0,
                'LambdaBl2'          : 4000,  #   223.503400817 obtaned from tov repo
                'domain'             : 0,      #Set 1 for FD. Default = 0
                'ecc'                : e,     # Eccentricity. Default = 0.
                'arg_out'            : "yes",      #Output hlm/hflm. Default = 0
                'use_mode_lm'        : k,      #List of modes to use/output through EOBRunPy
                'output_lm'          : k,      #List of modes to print on file
                'srate_interp'       : 4096.,  #srate at which to interpolate. Default = 4096.
                'use_geometric_units': "yes",   #output quantities in geometric units. Default = "yes"
                'initial_frequency'  : 0.002,    #in Hz if use_geometric_units = 0, else in geometric units
                'time_shift_TD'      : "yes",
                'interp_uniform_grid': "yes"   #interpolate mode by mode on a uniform grid. Default = "no" (no interpolation)
            }
teb, hp, hcm, hlm, dyn = EOBRun_module.EOBRunPy(pars)
nu = q_to_nu(q)
Ah22   = hlm['1'][0]*nu
Phih22 = hlm['1'][1]
heb = Ah22 * np.exp(-1j * (Phih22))
Momg22 = phi_dot(teb, heb)

# MomegaID = 2 * self.M * omega_NR
#print('Omega for NR ID = ', Momg22/(2*M))


## BBH
pars['LambdaAl2'] = 0
pars['LambdaBl2'] = 0
tebb, hp, hcm, hlm, dynb = EOBRun_module.EOBRunPy(pars)
Ah22b   = hlm[mode][0]*nu
Phih22b = hlm[mode][1]
hebb = Ah22b * np.exp(-1j * (Phih22b))


cecc = '#b2182b'
ceccp = '#2166ac'
fig = plt.figure(figsize=(13,6))
plt.plot(tebb,np.real(hebb),color='grey',lw=3,alpha=0.7,label='BBH')
plt.plot(teb-teb[np.argmax(np.abs(heb))],np.real(heb),color=cecc,alpha=0.7,label='BHNS')
plt.legend()
#plt.xlim([-400,400])
#plt.tight_layout()
#plt.subplots_adjust(bottom=0.15,wspace=0.08, hspace=0.08)
#plt.savefig("../../fig/ecc_eccprec_wvf.pdf")
#print("mass ratio = ", q, ", spin = ", spins[ind[0]],", lambda = ",lambdas[ind[1]])
plt.show()

# orbits
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(dyn['phi'], dyn['r'],color=cecc,alpha=0.7,label='BHNS')
ax.plot(dynb['phi'], dynb['r'],color='grey',lw=3,alpha=0.7,label='BBH')
plt.tight_layout()
plt.show()

# binding energy
pks = find_peaks(dyn['MOmega'], distance=150) [0]
nu  = pars['q']/(1 + pars['q'])**2
Eb = (dyn['E'] - 1.)/nu
Ebb = (dynb['E'] - 1.)/nu

# Plot
fig, ax = plt.subplots()
ax.plot(dynb['Pphi'], Ebb, color='grey',lw=3,alpha=0.7,label='BBH')
ax.plot(dyn['Pphi'], Eb,color=cecc,alpha=0.7,label='BHNS')
ax.scatter(dyn['Pphi'][pks], Eb[pks], color='k', label=r'Peaks of $\Omega_{\rm orb}$')

ax.set_xlabel(r'$p_\varphi$')
ax.set_ylabel(r'$\hat{E}_b$')
plt.legend()
plt.tight_layout()
plt.show()

# orbital frequency evolution
fig = plt.figure(figsize=(13,6))
plt.plot(dynb['t'],dynb['MOmega_orb'],color='grey',lw=3,alpha=0.7,label='BBH')
plt.plot(dyn['t'],dyn['MOmega_orb'],color=cecc,alpha=0.7,label='BHNS')
plt.legend()
plt.show()