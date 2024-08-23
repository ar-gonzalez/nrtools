import numpy as np
import cmath
import math
from scipy.signal import tukey
from watpy.wave.wave import rinf_str_to_float
import re
from scipy import integrate 
import csv

# Constants
Msun_sec = 4.925794970773135e-06

def float_to_latex_sci(f, precision=1):
    """
    Convert a float to LaTeX scientific notation.
    
    :param f: The float number to convert.
    :param precision: Number of decimal places.
    :return: A string representing the float in LaTeX scientific notation.
    """
    if f == 0:
        return f'0'
    
    exponent = int(np.floor(np.log10(abs(f))))
    mantissa = round(f / 10**exponent, precision)
    
    return f"{mantissa} \\times 10^{{{exponent}}}"

def csv_reader(filename):
    """
    Read a CSV file using csv.DictReader
    """
    data = [] 
    with open(filename) as f:
        reader = csv.DictReader(f, delimiter=',')
        for line in reader:
            data.append(line)
    return data

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

def get_kappa2t(lam,m1,m2):
    q = m1/m2
    nu = q_to_nu(q)
    lamtilde = (8/13)*(1+7*nu-31*nu*nu-np.sqrt(1-4*nu)*(1+9*nu-11*nu*nu))*lam
    return 3*lamtilde / 16    

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

    return Px_all, Py_all, Pz_all, P_all

##############################################################
# Compute hp, hc for SXS waveforms (ligic similar to C code)
##############################################################
from math import factorial as fact
from scipy import optimize
from scipy import interpolate
KMAX   = 14

def modes_to_k(modes):
    """
    Map multipolar (l,m) -> linear index k
    """
    return [int(x[0]*(x[0]-1)/2 + x[1]-2) for x in modes]

def k_to_ell(k):
    LINDEX = [\
    2,2,\
    3,3,3,\
    4,4,4,4,\
    5,5,5,5,5,\
    6,6,6,6,6,6,\
    7,7,7,7,7,7,7,\
    8,8,8,8,8,8,8,8]
    return LINDEX[k]

def k_to_emm(k):
    MINDEX = [\
    1,2,\
    1,2,3,\
    1,2,3,4,\
    1,2,3,4,5,\
    1,2,3,4,5,6,\
    1,2,3,4,5,6,7,\
    1,2,3,4,5,6,7,8];
    return MINDEX[k]   

def spinsphericalharm(s, l, m, phi, i):  
  c = pow(-1.,-s) * np.sqrt( (2.*l+1.)/(4.*np.pi) )
  dWigner = c * wigner_d_function(l,m,-s,i)
  rY = np.cos(m*phi) * dWigner
  iY = np.sin(m*phi) * dWigner
  return rY + 1j*iY

def wigner_d_function(l,m,s,i):
  costheta = np.cos(i*0.5)
  sintheta = np.sin(i*0.5)
  norm = np.sqrt( (fact(l+m) * fact(l-m) * fact(l+s) * fact(l-s)) )
  ki = max( 0  , m-s )
  kf = min( l+m, l-s )
  dWig = 0.
  for k in range(int(ki), int(kf)+1):
    div = 1.0/( fact(k) * fact(l+m-k) * fact(l-s-k) * fact(s-m+k) )
    dWig = dWig+div*( pow(-1.,k) * pow(costheta,2*l+m-s-2*k) * pow(sintheta,2*k+s-m) )
  
  return (norm * dWig)

def get_A_NR(dict, N=0):
    """
    Compute the aplitude of the NR waveform, skipping the first N points
    """
    r = np.sqrt(dict[N:,1]**2 +  dict[N:,2]**2)
    return r

def get_p_NR(dict, N=0):
    """
    Compute the phase of the NR waveform, skipping the first N points
    """
    r = -np.unwrap(np.arctan2(dict[N:,2], dict[N:,1]))
    return r

def InterpAmpPhase(A, p, t, tnew):
    """
    Interpolate Amplitude(t) and phase(t) on a new time array tnew
    """
    Af   = interpolate.interp1d(t, A, fill_value='extrapolate')
    Pf   = interpolate.interp1d(t, p, fill_value='extrapolate')
    Anew = Af(tnew)
    Pnew = Pf(tnew)
    return Anew, Pnew

def ComputeHpHc(d, d_neg, d0, i, prf, activemode):
    Ylm_d     = {}
    Ylm_d_neg = {}
    Ylm_d_0   = {}
    ellmax    = 2
    
    for k in range(0, KMAX):
        ell = k_to_ell(k)
        m   = k_to_emm(k)
        Ylm_d[k]    = spinsphericalharm(-2, ell, m, np.pi/2, i)
        Ylm_d_neg[k]= spinsphericalharm(-2, ell,-m, np.pi/2, i)
        
        ellmax = ell
    
    for l in range(2, ellmax+1):
        Ylm_d_0[l]  = spinsphericalharm(-2, l, 0, 0., i)

    sumr=0 
    sumi=0
    # Loop over modes
    for k in range(0, KMAX): 
        
        if (not activemode[k]): 
            continue
        Aki    = prf*d[k][0]
        cosPhi = np.cos(d[k][1])
        sinPhi = np.sin(d[k][1])
        sumr   = sumr+ Aki*(cosPhi*np.real(Ylm_d[k]) + sinPhi*np.imag(Ylm_d[k]))
        sumi   = sumi+ Aki*(cosPhi*np.imag(Ylm_d[k]) - sinPhi*np.real(Ylm_d[k])) 
          
	    # add m<0 modes
        Aki    = prf*d_neg[k][0]
        cosPhi = np.cos(d_neg[k][1])
        sinPhi = np.sin(d_neg[k][1])
        sumr   = sumr+ Aki*(cosPhi*np.real(Ylm_d_neg[k]) + sinPhi*np.imag(Ylm_d_neg[k]))
        sumi   = sumi+ Aki*(cosPhi*np.imag(Ylm_d_neg[k]) - sinPhi*np.real(Ylm_d_neg[k])) 
    
    # M = 0
    #for l in range(2, ellmax+1):
    #    if (d0[l] == None):
    #        continue
    #    Aki    = prf*d0[l][0]
    #    cosPhi = np.cos(d0[l][1])
    #    sinPhi = np.sin(d0[l][1])
    #    sumr   = sumr+ Aki*(cosPhi*np.real(Ylm_d_0[l]) + sinPhi*np.imag(Ylm_d_0[l]))
    #    sumi   = sumi+ Aki*(cosPhi*np.imag(Ylm_d_0[l]) - sinPhi*np.real(Ylm_d_0[l]))
      
    #h = h+ - i hx
    hp = sumr
    hc = -sumi
    
    return hp, hc

def ComputeHpHc_SXS(gw_sxs, lmax, i, prf, M, dT, N):
    d     = {}
    d_neg = {}
    d0    = {}
    Msuns  = 4.925491025543575903411922162094833998e-6
    # initialize d0 and actmodes:
    actmodes = np.zeros(KMAX)
    for l in range(2, k_to_ell(KMAX)+1):
        d0[l] = None
    
    # read all modes in dictionary
    for ell in range(2, lmax+1):
        for m in range(-ell, ell+1):
            if m == 0:
                continue #skip m = 0 modes

            ylm_str = "Y_l"+str(ell)+"_m"+str(m)+".dat" 
            gw_ext = gw_sxs["Extrapolated_N4.dir"][ylm_str]
            A      = get_A_NR(gw_ext, N)
            p      = get_p_NR(gw_ext, N)
            t_SI   = np.array(gw_ext[N:,0])*M*Msuns
            t_new  = np.arange(t_SI[0], t_SI[-1], dT)
            A, p   = InterpAmpPhase(A,  p,  t_SI, t_new)
            if m < 0:   #handle cases separately
                if m==-1:
                    continue
                m = -m
                k           = modes_to_k([[ell, m]])[0]
                d_neg[k]    = [A, p]
                actmodes[k] = 1

            elif m > 0:
                if m==1:
                    continue
                if m==2:
                    k           = modes_to_k([[ell, m]])[0]
                    d[k]        = [A, p]
                    actmodes[k] = 1
            else:
                d0[ell]     = [A, p]
    # from dictionary, compute hp,hc
    hp, hc = ComputeHpHc(d, d_neg, d0, i, prf, actmodes)
    return t_new, hp, hc

