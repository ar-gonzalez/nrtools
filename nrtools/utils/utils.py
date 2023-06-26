import numpy as np
import cmath
import math
from scipy.signal import tukey
from watpy.wave.wave import rinf_str_to_float
import re

# Constants
Msun_sec = 4.925794970773135e-06

## 3PN energy and ang. mom
# from https://journals.aps.org/prd/pdf/10.1103/PhysRevD.65.124009
# Eq. (3) and (4), for spins add Eq. (13)

def eb_3pn(mo,q):
    # mo is the angular orbital frequency: MOmega
    lam = -1987./3080. # not the tidal polarizability
    ost = 0 # omega_static in the paper
    nu = q_to_nu(q)
    term1 = 1 + (- 0.75 - nu/12)*pow(mo,2/3)
    term2 = (- 27/8 + nu*19/8 - nu*nu/24)*pow(mo,4/3)
    inside_term = 209323/4032 - np.pi*np.pi*205/96 - lam*110/9
    term3 = (- 675/64 + inside_term*nu - nu*nu*155/96 - nu*nu*nu*35/5184)*mo*mo
    return -0.5*nu*pow(mo,2/3)*(term1+term2+term3)

def j_3pn(mo,q):
    lam = -1987. / 3080.
    ost = 0 # omega_static in the paper    
    nu = q_to_nu(q)
    term1 = 1 + (1.5 + nu/6)*pow(mo,2/3)
    term2 = ( 27/8 - nu*19/8 + nu*nu/24)*pow(mo,4/3)
    inside_term = -209323/5040 + np.pi*np.pi*41/24 + lam*88/9
    term3 = ( 135/16 + inside_term*nu + nu*nu*31/24 + nu*nu*nu*7/1296)*mo*mo
    return nu*pow(mo,-1/3)*(term1+term2+term3)

## Misc
def get_id_gw_frequency_Hz(omega):
    # omega : orbital angular velocity from ID
    return omega / (math.pi * Msun_sec)

def get_id_gw_frequency_Momega22(omega, mtot):
    # omega : orbital angular velocity from ID
    # mtot : total mass of binary from ID
    return 2 * mtot * omega

def num_orbits_1pn(m1,m2,omega,mtot):
    # m1= BH_Christodoulou_mass_current
    # m2= NS_adm_mass
    # omega = angular velocity
    # mtot = ADM_mass
    nu = m1*m2/pow(m1+m2,2)
    return pow(mtot*omega,-5./3.)/(32.*nu)/(2.*np.pi)

def nu_to_q(nu):
    if nu==0:
        sol = 0.
    else: 
        a = nu
        b = 2*nu - 1
        c = nu

        d = b**2 - 4*a*c

        if d<0:
            sol_1 = (-b + cmath.sqrt(d))/(2*a)
            sol_2 = (-b - cmath.sqrt(d))/(2*a)
        else:
            sol_1 = (-b + math.sqrt(d))/(2*a)
            sol_2 = (-b - math.sqrt(d))/(2*a)
    
        if sol_1>0:
            sol = sol_1
        else:
            sol = sol_2

    return sol

def q_to_nu(q):
    return q/((1.+q)*(1.+q))

def get_rad(filename):
    name = re.match(r'R(\w+)_l(\d+)_m(\w+)_r(\w+).txt', filename)
    try:
        rad1   = rinf_str_to_float(name.group(4))
    except AttributeError:
        name2 = re.match(r'm(\d+)', name.group(4))
        rad1 = rinf_str_to_float(name2.group(1))
    return rad1


def windowing(h, alpha=0.1): 
   """ Perform windowing with Tukey window on a given strain (time-domain) 
       __________ 
       h    : strain to be tapered 
       alpha : Tukey filter slope parameter. Suggested value: alpha = 1/4/seglen 
       Only tapers beginning of wvf, to taper both, comment: window[len(h)//2:] = 1.
   """ 
   window = tukey(len(h), alpha) 
   #window[len(h)//2:] = 1. 
   wfact  = np.mean(window**2) 
   window = np.array(window) 
   return h*window, wfact
