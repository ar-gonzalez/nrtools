# Py utils
import os
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import argparse 
from scipy.signal import find_peaks
import sympy as sp

# teobresums run module 
import EOBRun_module

# watpy utils
from watpy.utils import num


# Obtain symmetric mass-ratio nu from q
def q_to_nu(q):
    return q/(1+q)**2

# Phase of the signal
def phase(signal):
    return -np.unwrap(np.angle(signal))

# omega of the signal
def phi_dot(time, signal):
    phase_to_diff = phase(signal)
    return num.diff1(time, phase_to_diff, pad=True)


class eob_qc_matter(object):

    def __init__(self, M, q, chi1, chi2, lam1, lam2):
        self.M = M
        self.q = q
        self.chi1 = chi1
        self.chi2 = chi2
        self.lam1 = lam1
        self.lam2 = lam2
        self.nu = q_to_nu(self.q)

    # EOB parameters
    def eob_pars(self):
        pars = {
        # System parametes, assuming aligned spins        
        'M'                  : self.M,      # Total mass
        'q'                  : self.q,      # Mass ratio m1/m2 > 1
        'chi1'               : self.chi1,     # Z component of chi_1
        'chi2'               : self.chi2,     # Z component of chi_2
        'LambdaAl2'          : 0.,     # Quadrupolar tidal parameter of body 1 (A)
        'LambdaBl2'          : self.lam2,     # Quadrupolar tidal parameter of body 2 (B)
        'ecc'                : 0.2,     # Eccentricity. Default = 0.
        #'ecc_freq'           : 1,      # Use periastron (0), average (1) or apastron (2) frequency forinitial condition computation. Default = 1

        # Initial conditions and output time grid
        'domain'             : 0,      # Time domain. EOBSPA is not available for eccentric waveforms!
        'initial_frequency'  : 0.001,     # in Hz if use_geometric_units = 0, else in geometric units
        'interp_uniform_grid': "yes",  # interpolate mode by mode on a uniform grid. Default = 0 (no interpolation)
        #'time_shift_TD'      : 0,      # t=0 at Merger (Default = 1)

        # Modes
        'use_mode_lm'        : [1],    # List of modes to use/output through EOBRunPy

        # Output parameters (Python)
        'arg_out'            : "yes",  # Output hlm/hflm. Default = 0
        }

        return pars

    # Run EOB
    def eob_run(self):
        self.t_eob, _, _, self.hlm_eob, self.dyn = EOBRun_module.EOBRunPy(self.eob_pars())


    # EOB useful quantities
    def eob_evolution_data(self):
        self.eob_run()
        self.A22_eob     = self.hlm_eob['1'][0]*self.nu
        self.p22_eob     = self.hlm_eob['1'][1]
        self.h22_eob     = self.A22_eob*np.exp(-1j*self.p22_eob)
        self.dot_p22_eob = phi_dot(self.t_eob, self.h22_eob)
        self.eob_idx     = np.argmax(self.A22_eob)
        self.t_eob       -= self.t_eob[self.eob_idx]

        self.eob_r       = self.dyn['r']
        self.eob_t_dyn   = self.dyn['t']
        self.eob_omg_dyn = self.dyn['MOmega']

        self.eob_data = {}
        self.eob_data['A22'] = self.A22_eob
        self.eob_data['p22'] = self.p22_eob
        self.eob_data['h22'] = self.h22_eob
        self.eob_data['Momg22'] = self.dot_p22_eob
        self.eob_data['r'] = self.eob_r * self.M
        self.eob_data['t'] = self.eob_t_dyn - np.abs(self.t_eob[0])
        self.eob_data['Momg'] = self.eob_omg_dyn

        return self.t_eob, self.eob_data


    # Plot
    def show_eob(self, omega):
        self.eob_evolution_data()

        MomegaID = 2 * self.M * omega
        t_ini_NR = np.interp(MomegaID, self.eob_data['Momg22'][:self.eob_idx], self.t_eob[:self.eob_idx])

        cut = (np.abs(self.t_eob - t_ini_NR)).argmin()
        peaks, _ = find_peaks(np.real(self.eob_data['h22'][cut:]), height=0)  # only peaks above 0
        peaks_neg, _ = find_peaks(-np.real(self.eob_data['h22'][cut:]), height=0)  # only peaks above 0

        r_ini = np.interp(t_ini_NR, self.eob_data['t'], self.eob_data['r'])
 
        if self.t_eob[cut:][peaks_neg[0]] < self.t_eob[cut:][peaks[0]]:
            orbits = peaks[0::2] 
        else:
            orbits = peaks[1::2]

        #orbits = peaks[0::2] 
        nr_orbits = len(orbits) - 1
        print()
        print("Number of full orbits:", nr_orbits)
        print("Initial separation guess from NR ID:", (self.M / omega**2)**(1/3))
        print("Initial EOB separation:", r_ini)
        print()

        
        fig, axs = plt.subplots(2, sharex = True, figsize = (8,6), gridspec_kw = {'height_ratios': [2,1]})
        fig.suptitle("EOB wave from initial SGRID Omega = {omega} || {nr_orbits} full orbits")

        axs[0].plot(self.t_eob, np.real(self.eob_data['h22']), c = '#d40000')
        axs[0].plot(self.t_eob, self.eob_data['A22'], c = '#d40000', lw = 0.5)
        for i, p in enumerate(orbits, start=0):
            axs[0].text(self.t_eob[cut:][p], np.real(self.eob_data['h22'][cut:])[p]+0.005, str(i), ha='center', va='bottom', fontsize=9, color='black')
        #axs[0].plot(self.t_eob[cut:][orbits], np.real(self.eob_data['h22'][cut:])[orbits], 'kx', label = f"Number of orbits: {nr_orbits}")
        axs[0].axvline(x = t_ini_NR, ls='--', c='k')
        axs[0].set_ylabel(r"$\Re h_{22}/M$, $A_{22}/M$")
        axs[0].set_ylim(-max(self.eob_data['A22'])*1.2, max(self.eob_data['A22'])*1.2)
        #axs[0].legend()

        axs[1].plot(self.t_eob, self.eob_data['Momg22'], c = '#d40000', label = 'GW')
        axs[1].plot(self.eob_data['t'], 2*self.eob_data['Momg'], 'k:', label = 'Orbital')
        axs[1].axvline(x = t_ini_NR, ls='--', c='k')
        axs[1].plot(t_ini_NR, MomegaID, 'gx', label = 'Initial NR $M\omega_{22}$')
        #axs[1].axhline(y = MomegaID, ls='--', c='b')
        axs[1].set_ylabel(r"$M\omega_{22}$")
        axs[1].set_xlabel(r"$u/M$")
        axs[1].set_xlim(t_ini_NR - 100, 100)
        axs[1].set_ylim(0.0, 0.25)
        axs[1].legend(loc = 'upper center')

        plt.tight_layout()
        plt.show()

        plt.plot(self.eob_data['t'], self.eob_data['r'])
        plt.axvline(x = t_ini_NR, ls='--', c='k')
        plt.plot(t_ini_NR, r_ini, 'gx', label = 'Guessed initial distance from EOB dynamics')
        plt.xlim(t_ini_NR - 100, 100)
        plt.ylim(0, 80)
        plt.xlabel("$u/M$")
        plt.ylabel("EOB separation")
        plt.legend()
        plt.tight_layout()
        plt.show()


if __name__ == "__main__": 
   
    parser = argparse.ArgumentParser(description='Script to check EOB waveform from greedy table params and SGRID Omega to estimate the number of orbits and eventually change the initial separation')

    parser.add_argument('-omg', '--Omega_ID', type=float, required=True, help='Omega par from SGRID') # 0.005

    args = parser.parse_args()

    q= 2.5
    m2 = 1.2498744738278933 #MS1b
    m1 = m2*q 
    lam1= 0 
    lam2=2312.470950148368
    chi1=0.75
    chi2 = 0

    M = m1 + m2

    eob = eob_qc_matter(M, q, chi1, chi2, lam1, lam2)
    eob.show_eob(args.Omega_ID)

    # estimate value of irreducible mass of the BH (needed as input for Elliptica)
    # m1 above is the christodolou mass of the BH

    # ---- define symbols ----
    x = sp.symbols('x')

    # ---- constants ----
    b = m1
    a = chi1

    # ---- define equation ----
    expr = x**4 - b**2 * x**2 + (1/4) * b**4 * a**2 # from definition of MChr

    # ---- solve equation ----
    solutions = sp.solve(expr, x)

    # ---- print results ----
    print(solutions)



