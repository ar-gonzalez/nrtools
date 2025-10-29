import numpy as np
from nrtools.initialdata.initialdata import *
from nrtools.utils.utils import *
import matplotlib.pyplot as plt
import EOBRun_module
from watpy.wave.gwutils import *
from watpy.utils.units import MSun_sec, MSun_meter
from watpy.wave.wave import wfile_get_detrad, wfile_parse_name
import matplotlib

matplotlib.rcParams['text.usetex']= True
matplotlib.rcParams['font.serif']= 'Palatino' 
matplotlib.rcParams['font.size']= 15 #28
import argparse
import os

def modes_to_k(modes):
    """
    Map multipolar (l,m) -> linear index k
    """
    return [int(x[0]*(x[0]-1)/2 + x[1]-2) for x in modes]

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description='Compares 1 simulation with its EOB waveform')

    # Add arguments
    parser.add_argument('-i','--input', help='ID name')
    parser.add_argument('-e','--evolution', help='Evolution name')
    parser.add_argument('-c','--cluster', help='Cluster where it is stored')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Access the values of the arguments
    simname = args.input
    evo = args.evolution
    cluster = args.cluster

    return simname, evo, cluster

if __name__ == '__main__':
    simname, evo, cluster = main()

    if cluster=='ARA':
        basedir = '/beegfs/mo63poy/BHNS_Elliptica'
        id_exe = '/home/mo63poy/BHNS_Elliptica/Elliptica/Exe/elliptica'
        ev_path = '/home/mo63poy/BNS_BAM/BAM'
    elif cluster=='DRACO':
        basedir = '/home/mo63poy/BHNS_INITIAL_DATA'
        id_exe = '/home/mo63poy/Elliptica/Exe/elliptica'
        ev_path = '/home/mo63poy/BAM'
    elif cluster=='home': # my PC
        basedir = '/home/agonzalez/Documents/PhD/NR/BHNS_evs/evolutions'
    else:
        print('Error: cluster name unknown')

    simpath = os.path.join(basedir,simname)
    bhns_id = Initial_Data(path=simpath,params=None, id_exe=id_exe)
    id_output = bhns_id.ou
    id_output.read_md()
    dic = id_output.id_dic
    chi1z = float(dic['BH_chi_z_current'])
    chi1y = float(dic['BH_chi_y_current'])
    chi1x = float(dic['BH_chi_x_current'])

    EOS = simname.split('_')[0]
    Mbh, Mg, M = id_output.get_msun_masses()

    q = Mbh/Mg
    nu = q_to_nu(q)

    lam2, _ = id_output.get_tidal_params()
    
    # Get NR Data
    ev1_path = os.path.join(simpath,evo)
    core_out = os.path.join(ev1_path,'CoReDB')

    #############################################################
    # Plotting all (l,m) modes available with (2,2) alignment
    #############################################################
    hlm_list = sorted([i for i in os.listdir(core_out) if i.startswith('Rh_l2_m')])
    #fig, ax = plt.subplots(4,1, figsize=(10,8))
    i = np.pi/3 # inclination
    phi = np.pi/2
    PC_SI  = 3.085677581491367e+16 # m
    MPC_SI = 1e6 * PC_SI

    umm = np.loadtxt(fname=os.path.join(core_out,hlm_list[0]), comments='#', usecols=(0), unpack=True)

    h = np.zeros_like( 1j* umm  )
    for file in hlm_list:
        vlmr = wfile_parse_name(file)
        var, l, m, rad, _ = vlmr
        print('m = ',m)
        uMlm, rhlm, irhlm, momg, aMlm, philm, tlm = np.loadtxt(fname=os.path.join(core_out,file), comments='#', usecols=(0,1,2,3,4,5,6), unpack=True)
        time = uMlm*MSun_sec()
        distance = rad*MPC_SI
        amplitude_prefactor = M * MSun_meter() / distance
        
        Alm = amplitude_prefactor*aMlm
        hlm = Alm * np.exp( - 1j * philm )

        ylm = spinsphericalharm(-2,l,m,phi,i)
        yl_m = spinsphericalharm(-2,l,-m,phi,i)
        
        h += hlm * ylm
        #plt.plot(uMlm,np.abs(h),label=str(m))
        #plt.legend()
        #plt.show()


    
    h0 = np.zeros_like( 1j* umm  )
    for file in hlm_list:
        vlmr = wfile_parse_name(file)
        var, l, m, rad, _ = vlmr
        if m==0:
            continue
        else:
            print('m = ',m)
            uMlm, rhlm, irhlm, momg, aMlm, philm, tlm = np.loadtxt(fname=os.path.join(core_out,file), comments='#', usecols=(0,1,2,3,4,5,6), unpack=True)
            time = uMlm*MSun_sec()
            distance = rad*MPC_SI
            amplitude_prefactor = M * MSun_meter() / distance

            Alm = amplitude_prefactor*aMlm
            hlm = Alm * np.exp( - 1j * philm )

            ylm = spinsphericalharm(-2,l,m,phi,i)
            yl_m = spinsphericalharm(-2,l,-m,phi,i)

            h0 += hlm * ylm

    hp = np.real(h)
    hc = -1* np.imag(h)
    al2 = np.abs(h)
    hp0 = np.real(h0)
    hc0 = -1* np.imag(h0)
    al0 = np.abs(h0)
    plt.plot(uMlm,al2,label=r'Amp $\ell=2$')
    #plt.plot(uMlm,hc,label='hcross')
    plt.plot(uMlm,al0,'--',label=r'Amp wo $m=0$')
    plt.legend()
    plt.show()

    plt.plot(uMlm,al2-al0)
    plt.yscale('log')
    plt.show()
        
