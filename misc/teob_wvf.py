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

def modes_to_k(modes):
    """
    Map multipolar (l,m) -> linear index k
    """
    return [int(x[0]*(x[0]-1)/2 + x[1]-2) for x in modes]


modes = [[2,1],[2,2],[3,2],[3,3],[4,4]]
k = modes_to_k(modes)
print(k)

mode = str(modes_to_k([[2,2]])[0])

q = 6
mns = 1.4
mbh = q*mns
M = mbh+mns
e = 0.2

pars = {
                'M'                  : M,
                'q'                  : q,
                'chi1x'              : 0,
                'chi1y'              : 0,
                'chi1z'              : 0.0,
                'chi2'               : 0.,
                'LambdaAl2'          : 0,
                'LambdaBl2'          : 500,  #   223.503400817 obtaned from tov repo
                'domain'             : 0,      #Set 1 for FD. Default = 0
                'ecc'                : 0,     # Eccentricity. Default = 0.
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
teb = teb - teb[np.argmax(Ah22)]

Ah33   = hlm['4'][0]*nu
Phih33 = hlm['4'][1]
heb33 = Ah33 * np.exp(-1j * (Phih33))


## BBH
pars['LambdaAl2'] = 0
pars['LambdaBl2'] = 0
tebb, hp, hcm, hlm, dynb = EOBRun_module.EOBRunPy(pars)
Ah22b   = hlm['1'][0]*nu
Phih22b = hlm['1'][1]
hebb = Ah22b * np.exp(-1j * (Phih22b))



plt.plot(tebb, np.real(hebb), label='BBH', color='gray')
plt.plot(teb, np.real(heb),label='BHNS 22')
plt.plot(teb, np.real(heb33),label='BHNS 33')
plt.xlim([10,teb[-1]])
plt.legend()
plt.show()

#plt.plot(teb,Ah33/Ah22)
plt.plot(teb,Ah22)
plt.plot(teb,Ah33)
plt.xlim([10,teb[-1]])
plt.show()