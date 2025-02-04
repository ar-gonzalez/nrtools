import numpy as np
import cmath
import math
from scipy.signal.windows import tukey
from watpy.wave.wave import rinf_str_to_float
import re
from scipy import integrate 
import csv
from scipy.linalg import eig, norm

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


def windowing(h, alpha=0.1,taper_beginning=False, taper_end=False): 
   """ Perform windowing with Tukey window on a given strain (time-domain) 
       __________ 
       h    : strain to be tapered 
       alpha : Tukey filter slope parameter. Suggested value: alpha = 1/4/seglen 
       if taper_beginning = taper_end = False it tapers both ends of wvf
   """ 
   window = tukey(len(h), alpha) 
   if taper_beginning:
       window[len(h)//2:] = 1.
   elif taper_end:
       window[:len(h)//2] = 1.
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

##########################################################
# PRECESSING WVFS TOOLS
# shamelessly taken from pyart
# ######################################################## 
def calc_initial_jframe(J,L, mwave):
    """
    Rotate multipoles wlm to an initial frame aligned with the
    total angular momentum J = L + S
    J: here the initial ADM angular momentum from the ID
    L: same as above but for the orb.ang.mom
    mwave: core mwaves object
    """
    #J = dyn['id']['J0']
    #L = dyn['id']['L0']

    # normalize J and get the rotation angles
    J_norm = norm(J)
    JJ     = J/J_norm
    thetaJ = np.arccos(JJ[2])
    phiJ   = np.arctan2(JJ[1],JJ[0])

    # euler angles for the rotation
    beta  = -thetaJ
    gamma = -phiJ

    # we have one additional dof to set. We choose it such that Lx = 0
    LL    = rotate3(L, 0, beta, gamma)
    psiL  = np.arctan2(LL.T[1], LL.T[0])
    alpha = -psiL

    euler_ang = np.array([alpha, beta, gamma])

    # Rotate the multipoles
    new_wvf = {}
    rad = mwave.radii[-1]
    for lm in mwave.modes:
        ell,emm = lm
        # find relevant multipoles to rotate
        same_ells = [ (l,m) for (l,m) in mwave.modes if l==ell ]
        #u = mwave.get(ell,emm,rad).time_ret()
        same_ells_dict = { (l,m): [mwave.get(l,m,rad).time_ret(), np.real(mwave.get(l,m,rad)), np.imag(mwave.get(l,m,rad))] 
                           for (l,m) in same_ells 
                        }
        # perform rotation
        rotd_hlm = rotate_wfarrs_at_all_times(ell, emm, same_ells_dict, euler_ang)
        new_wvf[(ell, emm)] = rotd_hlm
    
    return new_wvf

def rotate3_axis(vector,theta=0., axis = [0,0,1]):
    """
    Rotate a 3 vector around a provided axis of an angle theta
    """
    from scipy.spatial.transform import Rotation

    zaxis    = np.array(axis)
    r        = Rotation.from_rotvec(theta*zaxis)
    vector_r =r.apply(vector)
    return vector_r


# Given dictionary of multipoles all with the same l, calculate the roated multipole with (l,mp)
def rotate_wfarrs_at_all_times( l,                          # the l of the new multipole (everything should have the same l)
                                m,                          # the m of the new multipole
                                like_l_multipoles_dict,     # dictionary in the format { (l,m): array([domain_values,re,img]) }
                                euler_alpha_beta_gamma,
                                ref_orientation = None ):           

    '''
    Given dictionary of multipoles all with the same l, calculate the roated multipole with (l,mp).
    Key reference -- arxiv:1012:2879
    Based on LL,EZH 2018
    '''
    # unpack the euler angles
    alpha,beta,gamma = euler_alpha_beta_gamma

    # Handle the default behavior for the reference orientation
    if ref_orientation is None:
        ref_orientation = np.ones(3)

    # Apply the desired reflection for the reference orientation. 
    # NOTE that this is primarily useful for BAM run which have an atypical coordinate setup if Jz<0
        
    gamma *= np.sign( ref_orientation[-1] )
    alpha *= np.sign( ref_orientation[-1] )

    new_plus  = 0
    new_cross = 0

    for lm in like_l_multipoles_dict:
        # See eq A9 of arxiv:1012:2879
        l,mp = lm
        old_wfarr = like_l_multipoles_dict[lm]

        d   = wdelement(l,m,mp,alpha,beta,gamma)
        a,b = d.real,d.imag
        p   = old_wfarr[1]
        c   = old_wfarr[2]

        new_plus  += a*p - b*c
        new_cross += b*p + a*c

    # Construct the new waveform array

    return  {'real': new_plus, 
            'imag' : new_cross, 
            'A'    : np.sqrt(new_plus**2+new_cross**2),
            'p'    : np.arctan2(new_cross,new_plus) }

# Rotate a 3 vector using Euler angles
def rotate3(vector,alpha,beta,gamma,invert=False):
    '''
    Rotate a 3 vector using Euler angles under conventions defined at:
    https://en.wikipedia.org/wiki/Euler_angles
    https://en.wikipedia.org/wiki/Rotation_matrix

    Science reference: https://arxiv.org/pdf/1110.2965.pdf (Appendix)

    Specifically, the Z1,Y2,Z3 ordering is used: https://wikimedia.org/api/rest_v1/media/math/render/svg/547e522037de6467d948ecf3f7409975fe849d07

    *  alpha represents a rotation around the z axis
    *  beta represents a rotation around the x' axis
    *  gamma represents a rotation around the z'' axis

    NOTE that in order to perform the inverse rotation, it is *not* enough to input different rotation angles. One must use the invert=True keyword. 
    This takes the same angle inputs as the forward rotation, but correctly applies the transposed rotation matricies in the reversed order.

    spxll'18
    '''

    # Import usefuls
    from numpy import cos,sin,array,dot,ndarray,vstack

    # Hangle angles as arrays
    angles_are_arrays = isinstance(alpha,np.ndarray) and isinstance(beta,np.ndarray) and isinstance(gamma,np.ndarray)
    if angles_are_arrays:
        # Check for consistent array shapes
        if not ( alpha.shape == beta.shape == gamma.shape ):
            # Let the people know and halt
            print( 'input angles as arrays must have identical array shapes' )

    # Validate input(s)
    if isinstance(vector,(list,tuple,ndarray)):
        vector = array(vector)
    else:
        print('first input must be iterable compatible 3D vector; please check')


    # Rotation around z''
    Ra = array( [
                    [cos(alpha),-sin(alpha),0],
                    [sin(alpha),cos(alpha),0],
                    [0,0,1]
        ] )

    # Rotation around y
    Rb = array( [
                    [cos(beta),0,sin(beta)],
                    [0,1,0],
                    [-sin(beta),0,cos(beta)]
        ] )

    # Rotation around z
    Rg = array( [
                    [cos(gamma),-sin(gamma),0],
                    [sin(gamma),cos(gamma),0],
                    [0,0,1]
        ] )

    # Perform the rotation
    # ans = (  Ra * ( Rb * ( Rg * vector ) )  )
    # NOTE that this is the same convention of equation A9 of Boyle et al : https://arxiv.org/pdf/1110.2965.pdf
    R = dot(  Ra, dot(Rb,Rg)  )
    if invert: R = R.T
    ans = dot( R, vector )

    # If angles are arrays, then format the input such that rows in ans correspond to rows in alpha, beta and gamma
    if angles_are_arrays:
        ans = vstack( ans ).T

    return ans


# Given a mwave, calculate the Euler angles corresponding to a co-precessing frame
# Taken from nrutils_dev, credits to Lionel London
def calc_coprecessing_angles(mwave,        # mwave object from watpy
                             domain_vals=None,      # The time or freq series for multipole data
                             ref_orientation=None,  # e.g. initial J; used for breaking degeneracies in calculation
                             return_xyz=False,
                             safe_domain_range=None):
    '''
    Key referece: https://arxiv.org/pdf/1304.3176.pdf
    Secondary ref: https://arxiv.org/pdf/1205.2287.pdf
    INPUT
    ---
    multipole_dict,       # dict of multipoles { ... l,m:data_lm ... }
    t,                    # The time series corresponding to multipole data; needed
                            only to calculate gamma; Optional
    OUTPUT
    ---
    alpha,beta,gamma euler angles as defined in https://arxiv.org/pdf/1205.2287.pdf
    AUTHOR
    ---
    Lionel London (spxll) 2017
    '''
    from scipy.linalg import eig, norm
    
    # Handle optional input
    if ref_orientation is None:
        ref_orientation = np.ones(3)
        
    #-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-#
    # Enforce that multipole data is array typed with a well defined length
    #-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-#
    mw = mwave
    rad = mw.radii[-1]
    y = {}
    h = {}
    for lm in mw.modes:
        l = lm[0]
        m = lm[1]
        w = mw.get(l=l,m=m,r=rad)
        h[lm] = w.h
        y[l, m] = np.array(h[lm])

    #-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-#
    # Calculate the emission tensor corresponding to the input data
    #-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-#
    L = calc_Lab_tensor(mw)

    #-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-#
    # Compute the eigenvectors and values of this tensor
    #-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-#

    # NOTE that members of L have the same length as each y[l,m]; the latter has been
    # forced to always have a length above

    # Initialize idents for angles. 
    # NOTE that gamma will be handled below
    alpha, beta = [], []
    X, Y, Z = [], [], []
    old_dom_dex = None

    # For all multipole instances
    ref_x, ref_y, ref_z = None, None, None
    flip_z_convention = False
    for k in range(len(L[0, 0, :])):

        # Select the emission matrix for this instance, k
        _L = L[:, :, k]

        # Compute the eigen vals and vecs for this instance
        vals, vec = eig(_L)

        # Find the dominant direction's index
        dominant_dex = np.argmax(vals)
        if old_dom_dex is None:
            old_dom_dex = dominant_dex
        if old_dom_dex != dominant_dex:
            # print dominant_dex
            old_dom_dex = dominant_dex

        # Select the corresponding vector
        dominant_vec = vec[:, dominant_dex]

        # There is a z axis degeneracy that we will break here
        # by imposing that the z component is always consistent with the initial L
        if not flip_z_convention:
            if np.sign(dominant_vec[-1]) == -np.sign(ref_orientation[-1]):
                dominant_vec *= -1
        else:
            if np.sign(dominant_vec[-1]) ==  np.sign(ref_orientation[-1]):
                dominant_vec *= -1

        # dominant_vec *= sign(domain_vals[k])

        # Extract the components of the dominant eigenvector
        _x, _y, _z = dominant_vec

        # Store reference values if they are None
        if ref_x == None:
            ref_x = _x
            ref_y = _y
            ref_z = _z
        else:
            if (ref_x*_x < 0) and (ref_y*_y < 0):
                _x *= -1
                _y *= -1
                _x *= -1

        # Store unit components for reference in the next iternation
        ref_x = _x
        ref_y = _y
        ref_z = _z

        # Look for and handle trivial cases
        if abs(_x)+abs(_y) < 1e-8:
            _x = _y = 0

        X.append(_x)
        Y.append(_y)
        Z.append(_z)

    # Look for point reflection in X
    X = reflect_unwrap(np.array(X))
    Y = np.array(Y)
    Z = np.array(Z)

    # 3-point vector reflect unwrapping
    # print safe_domain_range
    tol = 0.1
    if safe_domain_range is None:
        tol2 = 1e-20
        safe_domain_range = [min(abs(domain_vals)) - tol2, max(abs(domain_vals)) + tol2] 
    safe_domain_range = np.array(safe_domain_range)
    for k in range(len(X))[1:-1]:
        if k > 0 and k < (len(domain_vals)-1):

            if (abs(domain_vals[k]) > min(abs(safe_domain_range))) and (abs(domain_vals[k]) < max(abs(safe_domain_range))):

                left_x_has_reflected = abs(X[k]+X[k-1]) < tol*abs(X[k-1])
                left_y_has_reflected = abs(Y[k]+Y[k-1]) < tol*abs(X[k-1])

                right_x_has_reflected = abs(X[k]+X[k+1]) < tol*abs(X[k])
                right_y_has_reflected = abs(Y[k]+Y[k+1]) < tol*abs(X[k])

                x_has_reflected = right_x_has_reflected or left_x_has_reflected
                y_has_reflected = left_y_has_reflected or right_y_has_reflected

                if x_has_reflected and y_has_reflected:
                    
                    if left_x_has_reflected:
                        X[k:] *= -1
                    if right_x_has_reflected:
                        X[k+1:] *= -1

                    if left_y_has_reflected:
                        Y[k:] *= -1
                    if right_y_has_reflected:
                        Y[k+1:] *= -1

                    Z[k:] *= -1

    # Make sure that imag parts are gone
    X = np.double(X)
    Y = np.double(Y)
    Z = np.double(Z)

    #################################################
    # Reflect Y according to nrutils conventions    #
    Y *= -1                                         #
    #################################################

    a = np.array(ref_orientation)/norm(ref_orientation)
    B = np.array([X, Y, Z]).T
    b = (B.T/norm(B, axis=1))
    
    # Here we define a test quantity that is always sensitive to each dimension.
    # NOTE that a simple dot product does not have this property if eg 
    # a component of the reference orientation is zero. 
    # There is likely a better solution here.
    test_quantity = sum( [ a[k]*b[k] if a[k] else b[k] for k in range(3) ] )

    mask = (domain_vals >= min(safe_domain_range)) & (
        domain_vals <= max(safe_domain_range))
    if 1*(test_quantity[mask][0]) < 0:
        print('flipping manually for negative domain')
        X = -X
        Y = -Y
        Z = -Z

    # Calculate Angles

    alpha = np.arctan2(Y, X)
    beta  = np.arccos(Z)

    # Make sure that angles are unwrapped
    alpha = np.unwrap(alpha)
    beta  = np.unwrap(beta)

    # Calculate gamma (Eq. A4 of of arxiv:1304.3176)
    if len(alpha) > 1:
        k = 1
        gamma = - spline_antidiff(domain_vals, np.cos(beta)
                                  * spline_diff(domain_vals, alpha, k=k), k=k)
        gamma = np.unwrap(gamma)
        # Enforce like integration constant for neg and positive frequency gamma; 
        # this assumes time series will not have negative values (i.e. the code should work for TD and FD cases)
        neg_mask = domain_vals < 0
        _mask = (-domain_vals) > 0.01
        mask_ = domain_vals > 0.01
        if sum(neg_mask):
            gamma[neg_mask] = gamma[neg_mask] - \
                gamma[_mask][-1] + gamma[mask_][0]
    else:
        # NOTE that this is the same as above, but here we're choosing an integration constant such that 
        # the value is zero. Above, no explicit integration constant is chosen.
        gamma = 0

    # Return answer
    if return_xyz == 'all':
        return alpha, beta, gamma, X, Y, Z
    elif return_xyz:
        return X, Y, Z
    else:
        return alpha, beta, gamma

# Calculate the emission tensor given a dictionary of multipole data
# Taken from nrutils_dev, credits to Lionel London
def calc_Lab_tensor(mwave):
    '''
    Given a dictionary of multipole moments (single values or time series)
    determine the emission tensor, <L(aLb)>.

    See:
    - https://arxiv.org/pdf/1304.3176.pdf
    - https://arxiv.org/pdf/1205.2287.pdf
    '''

    # Rename multipole_dict for short-hand
    mw = mwave
    rad = mw.radii[-1]
    x = {}
    y = {}
    h = {}
    for lm in mw.modes:
        l = lm[0]
        m = lm[1]
        w = mw.get(l=l,m=m,r=rad)
        h[lm] = w.h
        x[l, m] = np.array(h[lm])
        x[l,m,'conj'] = x[l,m].conj()

    y = x
    # Check type of dictionary values and pre-allocate output
#    if isinstance( y[2,2], (float,int,complex) ):
#        L = np.zeros( (3,3), dtype=complex )
#    elif isinstance( y[2,2], np.ndarray ):
    L = np.zeros( (3,3, len(y[2,2])), dtype=complex )
#    else:
#        print('Dictionary values of handled type; must be float or array')

    # define lambda function for useful coeffs
    c = lambda l,m: np.sqrt( l*(l+1) - m*(m+1) ) if abs(m)<=l else 0

    # Compute tensor elements (Eqs. A1-A2 of https://arxiv.org/pdf/1304.3176.pdf)
    I0,I1,I2,Izz = np.zeros_like(y[2,2]), np.zeros_like(y[2,2]), np.zeros_like(y[2,2]), np.zeros_like(y[2,2])

    # Sum contributions from input multipoles
    for lm in mw.modes:
        l = lm[0]
        m = lm[1]
        # Eq. A2c
        I0 += 0.5 * ( l*(l+1)-m*m ) * y[l,m] * y[l,m,'conj']
        # Eq. A2b
        I1 += c(l,m) * (m+0.5) * ( y[l,m+1,'conj'] if (l,m+1) in y else 0 ) * y[l,m]
        # Eq. A2a
        I2 += 0.5 * c(l,m) * c(l,m+1) * y[l,m] * ( y[l,m+2,'conj'] if (l,m+2) in y else 0 )
        # Eq. A2d
        Izz += m*m * y[l,m] * y[l,m,'conj']

    # Compute the net power (amplitude squared) of the multipoles
    N = sum( [ y[l,m] * y[l,m,'conj'] for l,m in mw.modes ] ).real

    # Populate the emission tensor ( Eq. A2e )
    # Populate asymmetric elements
    L[0,0] = I0 + I2.real
    L[0,1] = I2.imag
    L[0,2] = I1.real
    L[1,1] = I0 - I2.real
    L[1,2] = I1.imag
    L[2,2] = Izz
    # Populate symmetric elements
    L[1,0] = L[0,1]
    L[2,0] = L[0,2]
    L[2,1] = L[1,2]

    # Normalize
    N[ N==0 ] = min( N[N>0] )
    L = L.real / N

    return L

def reflect_unwrap( vec ):
    '''
    Reflect points in an array
    '''

    ans = np.array(vec)
    for k in range(len(vec)):
        if (k>0) and ( (k+1) < len(vec) ):
            l = vec[k-1]
            c = vec[k]
            r = vec[k+1]
            reflect = (np.sign(l)==np.sign(r)) and (np.sign(l)==-np.sign(c))
            if reflect:
                ans[k] *= -1
    return ans

def spline_antidiff(t,y,k=3,n=1):
    """
    Wrapper for InterpolatedUnivariateSpline antiderivative function
    """

    #
    from scipy.interpolate import InterpolatedUnivariateSpline as spline

    # Calculate the desired number of integrals
    ans = spline(t,y.real,k=k).antiderivative(n=n)(t) + ( 1j*spline(t,y.imag,k=k).antiderivative(n=n)(t) if isinstance(y[0],complex) else 0 )

    # Return the answer
    return ans

def spline_diff(t,y,k=3,n=1):
    """
    Wrapper for InterpolatedUnivariateSpline derivative function
    """

    #
    from numpy import sum
    from scipy.interpolate import InterpolatedUnivariateSpline as spline

    # Calculate the desired number of derivatives
    ans = spline(t,y.real,k=k).derivative(n=n)(t) \
          + ( 1j*spline(t,y.imag,k=k).derivative(n=n)(t) if (sum(abs(y.imag))!=0) else 0 )

    return ans

# Calculate Widger D-Matrix Element
def wdelement( ll,         # polar index (eigenvalue) of multipole to be rotated (set of m's for single ll )
               mp,         # member of {all em for |em|<=l} -- potential projection spaceof m
               mm,         # member of {all em for |em|<=l} -- the starting space of m
               alpha,      # -.
               beta,       #  |- Euler angles for rotation
               gamma ):    # -'

    #** James Healy 6/18/2012
    #** wignerDelement
    #*  calculates an element of the wignerD matrix
    # Modified by llondon6 in 2012 and 2014
    # Converted to python by spxll 2016
    #
    # This implementation apparently uses the formula given in:
    # https://en.wikipedia.org/wiki/Wigner_D-matrix
    #
    # Specifically, this the formula located here: 
    # https://wikimedia.org/api/rest_v1/media/math/render/svg/53fd7befce1972763f7f53f5bcf4dd158c324b55

    #
    if ( (type(alpha) is np.ndarray) and (type(beta) is np.ndarray) and (type(gamma) is np.ndarray) ):
        alpha,beta,gamma = alpha.astype(float), beta.astype(float), gamma.astype(float)
    else:
        alpha,beta,gamma = float(alpha),float(beta),float(gamma)

    coefficient = np.sqrt( fact(ll+mp)*fact(ll-mp)*fact(ll+mm)*fact(ll-mm))*np.exp( 1j*(mp*alpha+mm*gamma) )

    total = 0
    # find smin
    if (mm-mp) >= 0      :  smin = mm - mp
    else                 :  smin = 0
    # find smax
    if (ll+mm) > (ll-mp) : smax = ll-mp
    else                 : smax = ll+mm

    if smin <= smax:
        for ss in range(smin,smax+1):
            A = (-1)**(mp-mm+ss)
            A *= np.cos(beta/2)**(2*ll+mm-mp-2*ss)  *  np.sin(beta/2)**(mp-mm+2*ss)
            B = fact(ll+mm-ss) * fact(ss) * fact(mp-mm+ss) * fact(ll-mp-ss)
            total += A/B

    element = coefficient*total
    return element

##############################################################
# Compute hp, hc for SXS waveforms (logic similar to C code)
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

