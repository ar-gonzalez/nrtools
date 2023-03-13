import numpy as np
import cmath
import math

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
    term1 = 1 + (- 1.5 - nu/6)*pow(mo,2/3)
    term2 = ( 27/8 - nu*19/8 + nu*nu/24)*pow(mo,4/3)
    inside_term = -209323/5040 + np.pi*np.pi*41/24 + lam*88/9
    term3 = ( 135/16 + inside_term*nu + nu*nu*31/24 + nu*nu*nu*7/1296)*mo*mo
    return nu*pow(mo,-1/3)*(term1+term2+term3)

## Misc
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
