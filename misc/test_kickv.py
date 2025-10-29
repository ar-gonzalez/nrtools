import argparse
import os
from nrtools.initialdata.initialdata import *
from nrtools.evolution.evolution import *
from nrtools.utils.utils import get_rad
from watpy.wave.wave import wfile_get_mass, wfile_get_detrad
from watpy.wave.gwutils import *
from watpy.utils.units import *
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

matplotlib.rcParams['text.usetex']= True
matplotlib.rcParams['font.serif']= 'Palatino'
matplotlib.rcParams['font.size']= 15 #28

def get_stuff(bhns_ev):
    h22_file = sorted([i for i in os.listdir(bhns_ev.core_out) if i.startswith('Rh_l2_m2_r')])[-1]
    rad = int(wfile_get_detrad(os.path.join(bhns_ev.core_out,h22_file)))
    time22, amp22 = np.loadtxt(fname=os.path.join(bhns_ev.core_out,h22_file), comments='#', usecols=(6,4), unpack=True)
    tmrg = time22[np.argmax(amp22)]

    with open(os.path.join(bhns_ev.ou.out_inv_dir,pxfile+r+'.l'+lv), 'r') as file:
        content = file.read()
        #if 'nan' in content:
        bhns_ev.get_lin_momentum()
        try:
            Px, Py, Pz, P, time = np.loadtxt(fname=os.path.join(bhns_ev.core_out,'P_r00'+str(rad)+'.txt'), comments='#', usecols=(0,1,2,3,5), unpack=True)
        except:
            Px, Py, Pz, P, time = np.loadtxt(fname=os.path.join(bhns_ev.core_out,'P_r0'+str(rad)+'.txt'), comments='#', usecols=(0,1,2,3,5), unpack=True)
    return Px, Py, Pz, P, time

lvl = {
    6: '#ca0020',
    8: '#f4a582',
    10: '#92c5de',
    11: '#0571b0'
}

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description='compare different refinement levels of a BAM simulation')

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
    _, _, M = id_output.get_msun_masses()

    r = '1'
    lv = '0'

    c = 299792458.0
    pxfile = 'Px_r'
    pyfile = 'Py_r'

    resolution, lmax, flux = ev_folder_parser(evo)
    bhns_ev = Evolution(path=os.path.join(simpath,evo), ev_path=ev_path, initial_data=bhns_id, resolution=int(resolution), lmax=lmax, lmax2=6,flux=flux)
   
    # do things for single case
    bhns_ev.get_core_data()
    h22_file = sorted([i for i in os.listdir(bhns_ev.core_out) if i.startswith('Rh_l2_m2_r')])[-1]
    time22, amp22 = np.loadtxt(fname=os.path.join(bhns_ev.core_out,h22_file), comments='#', usecols=(6,4), unpack=True)
    tmrg = time22[np.argmax(amp22)]
    Px, Py, Pz, P, time = get_stuff(bhns_ev)
    #vel = (Px - 1j*Py)/M
    #v = np.abs(vel)
    v = P/M
    #v = np.sqrt(Px*Px + Py*Py + Pz*Pz) / M
    v = v*c/1000

    # total velocity
    #plt.plot(time,Px,label='Px') # color=lvl[11])#,label='Px')
    #plt.plot(time,Py,label='Py')
    #plt.plot(time,Pz,label='Pz')
    #plt.plot(time,P,'--',label='P')
    plt.plot(time,v)
    #plt.yscale('log')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.xlim([0,time[-1]])
    plt.xlabel(r't/M')
    plt.ylabel(r'$v^{\rm GW}_{\rm kick}~[km/s]$ ')
    plt.tight_layout()
    plt.show()
