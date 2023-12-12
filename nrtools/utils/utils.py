import numpy as np
import cmath
import math
from scipy.signal import tukey
from watpy.wave.wave import rinf_str_to_float
import re
from scipy import integrate 

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

# multipolar coefficients to get energetics from wvf
def mc_f(l,m):
    return np.sqrt(l*(l+1) - m*(m+1))
def mc_a(l,m):
    return np.sqrt((l-m)*(l+m+1))/(l*(l+1))
def mc_b(l,m):
    return np.sqrt(((l-2)*(l+2)*(l+m)*(l+m-1))/((2*l-1)*(2*l+1)))/(2*l)
def mc_c(l,m):
    return 2*m/(l*(l+1))
def mc_d(l,m):
    return np.sqrt(((l-2)*(l+2)*(l-m)*(l+m))/((2*l-1)*(2*l+1)))/l


def lin_momentum_from_wvf(h, doth, t, u, lmmodes):
    # From https://arxiv.org/pdf/0912.1285.pdf
    #* h[(l,m)]     : multipolar strain 
    #* doth[(l,m)]  : time-derivative of multipolar strain
    #* t            : time array
    #* lmmodes      : (l,m) indexes

    oo16pi  = 1./(16*np.pi)

    lmodes = [lm[0] for lm in lmmodes]
    mmodes = [lm[1] for lm in lmmodes]
    lmin = min(lmodes)
    if lmin < 2:
        raise ValueError("l>2")
    if lmin != 2:
        print("Warning: lmin > 2")
        
    mnfactor = np.ones_like(mmodes)

    dotP, P = {}, {}
    dotPz, dotPy, dotPx = {}, {}, {}
    Pz, Py, Px = {}, {}, {}

    dotP_all, P_all = np.zeros_like(t), np.zeros_like(t)
    dotPz_all, dotPy_all, dotPx_all = np.zeros_like(t), np.zeros_like(t), np.zeros_like(t)
    Pz_all, Py_all, Px_all = np.zeros_like(t), np.zeros_like(t), np.zeros_like(t)

    for k, (l,m) in enumerate(lmmodes):
        fact = mnfactor[k] * oo16pi
        dothlm1   = doth[(l,m+1)]   if (l,m+1)   in doth else 0*h[(l,m)]
        dothl_1m1 = doth[(l-1,m+1)] if (l-1,m+1) in doth else 0*h[(l,m)]
        dothl1m1  = doth[(l+1,m+1)] if (l+1,m+1) in doth else 0*h[(l,m)]
        dotl_1m   = doth[(l-1,m)]   if (l-1,m)   in doth else 0*h[(l,m)]
        dothl1m   = doth[(l+1,m)]   if (l+1,m)   in doth else 0*h[(l,m)]
        
        dotPxiy = 2.0 * fact * doth[(l,m)] * \
            (mc_a(l,m) * np.conj(dothlm1) + mc_b(l,-m) * np.conj(dothl_1m1) - mc_b(l+1,m+1) * np.conj(dothl1m1))
        dotPy[(l,m)] = np.imag(dotPxiy)
        dotPx[(l,m)] = np.real(dotPxiy)
        dotPz[(l,m)] = fact * np.imag( doth[(l,m)] * \
            (mc_c(l,m) * np.conj(doth[(l,m)]) + mc_d(l,m) * np.conj(dotl_1m) + mc_d(l+1,m) * np.conj(dothl1m)) )

        dotP[(l,m)] = np.sqrt(dotPx[(l,m)]**2 + dotPy[(l,m)]**2 + dotPz[(l,m)]**2)

        Pz[(l,m)] = integrate.cumtrapz(dotPz[(l,m)],t,initial=0)
        Py[(l,m)] = integrate.cumtrapz(dotPy[(l,m)],t,initial=0)
        Px[(l,m)] = integrate.cumtrapz(dotPx[(l,m)],t,initial=0)
        P[(l,m)]  = integrate.cumtrapz(dotP[(l,m)],t,initial=0)

        dotPz_all += dotPz[(l,m)]
        dotPy_all += dotPy[(l,m)]
        dotPx_all += dotPx[(l,m)]
        dotP_all  += dotP[(l,m)]

        Pz_all += Pz[(l,m)]
        Px_all += Px[(l,m)]
        Py_all += Py[(l,m)]
        P_all  += P[(l,m)]

    return Px, Py, Pz, P