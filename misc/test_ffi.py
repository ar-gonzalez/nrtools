###
# This script generates the CoRe wvf output
# and shows a plot of the strain
###
from nrtools.initialdata.initialdata import *
from nrtools.evolution.evolution import *
from watpy.wave.gwutils import fixed_freq_int_2
import os
import numpy as np
import argparse
import matplotlib.pyplot as plt
#import matplotlib

#matplotlib.rcParams['text.usetex']= True
#matplotlib.rcParams['font.serif']= 'Palatino' 
#matplotlib.rcParams['font.size']= 15 #28

## Read args
parser = argparse.ArgumentParser(description='Get Waveform from evolution, e.g. `python get_wvf.py -i id_name -e ev_name -c cluster_name`, either in ARA, DRACO, or HAWK')
# add arguments
parser.add_argument('-i', '--input', type=str, required=True, help='simulation name')
parser.add_argument('-e', '--evolution', type=str, required=True, help='evolution name')
parser.add_argument('-c', '--cluster', type=str, required=True, help='cluster')

# parse the arguments
args = parser.parse_args()

# access the arguments
simname = args.input
evo = args.evolution
cluster = args.cluster

if cluster=='ARA':
    basedir = '/beegfs/mo63poy/BHNS_Elliptica'
    id_exe = '/home/mo63poy/BHNS_Elliptica/Elliptica/Exe/elliptica'
    ev_path = '/home/mo63poy/BNS_BAM/BAM'
elif cluster=='DRACO':
    basedir = '/home/mo63poy/BHNS_INITIAL_DATA'
    id_exe = '/home/mo63poy/Elliptica/Exe/elliptica'
    ev_path = '/home/mo63poy/BAM'
elif cluster=='HAWK': # modify HERE!
        batchsys = 'pbs'
        basedir = '/zhome/academic/HLRS/xuj/xujapigo/BHNS_INITIAL_DATA'
        id_exe = '/zhome/academic/HLRS/xuj/xujapigo/Elliptica/Exe/elliptica'
        bam_path = '/zhome/academic/HLRS/xuj/xujapigo/BAM'
        workspace = '/lustre/hpe/ws10/ws10.3/ws/xujapigo-bhns'

        evo_path = os.path.join(workspace,simname) # where the simulation is running (inside our hawk workspace)
else:
    print('Error: cluster name unknown')

simpath = os.path.join(basedir,simname)
if cluster=='HAWK':
    evoo_path = os.path.join(evo_path,evo) # path to the workspace where the simulation is
else:
    evoo_path = os.path.join(simpath,evo)

## Do evolution(s)
# Initialize ID object
bhns_id = Initial_Data(path=simpath,params=None, id_exe=id_exe)

#ev_folders = [i for i in os.listdir(simpath) if os.path.isdir(os.path.join(simpath,i)) and i.startswith('bam')]
lmax2 = 6
print('\n>> ', evo)
resolution, lmax, flux = ev_folder_parser(evo)
bhns_ev = Evolution(path=evoo_path, ev_path=ev_path, initial_data=bhns_id, resolution=int(resolution), lmax=int(lmax), lmax2=lmax2,flux=flux) 


ev_output = bhns_ev.ou
id_output = bhns_id.ou
    

## cordb
files = os.listdir(bhns_ev.core_out)
filename = [i for i in files if i.startswith('Rh_l2_m1_r')][0]
time44, rh44 = np.loadtxt(fname=os.path.join(bhns_ev.core_out,filename), comments='#', usecols=(0,1), unpack=True)

## individual
mbh, mns, mtot = id_output.get_msun_masses()
_, Momg22 = id_output.get_gw_freqs()
wm = ev_output.get_mp_Rpsi4(mtot,Momg22) #bhns_ev.get_core_wm_object() #ev_output.get_mp_Rpsi4(mtot,Momg22)

h44 = wm.get(l=2,m=1)
m = 1
omega = 0.0008
fcut = 2 * omega / m 
h44.h = h44.get_strain(fcut=fcut, win=1.)

#h44.write_to_txt('h', bhns_ev.core_out)


##### print EOB prediction too
import EOBRun_module
from watpy.wave.gwutils import *

def modes_to_k(modes):
    """
    Map multipolar (l,m) -> linear index k
    """
    return [int(x[0]*(x[0]-1)/2 + x[1]-2) for x in modes]

modes = [[2,1],[2,2],[3,2],[3,3],[4,4]]
k = modes_to_k(modes)
    
id_output.read_md()
dic = id_output.id_dic
chi1z = float(dic['BH_chi_z_current'])
chi1y = float(dic['BH_chi_y_current'])
chi1x = float(dic['BH_chi_x_current'])

Mbh, Mg, M = id_output.get_msun_masses()

q = Mbh/Mg
nu = q_to_nu(q)

lam2, _ = id_output.get_tidal_params()

pars = {
                'M'                  : M, 
                'q'                  : q,
                'chi1x'              : chi1x,
                'chi1y'              : chi1y,
                'chi1z'              : chi1z,
                'chi2'               : 0.,
                'LambdaAl2'          : 0.,
                'LambdaBl2'          : lam2,  #   223.503400817 obtaned from tov repo
                'domain'             : 0,      #Set 1 for FD. Default = 0
                'arg_out'            : "yes",      #Output hlm/hflm. Default = 0
                'use_mode_lm'        : k,      #List of modes to use/output through EOBRunPy
                'output_lm'          : k,      #List of modes to print on file
                'srate_interp'       : 4096.,  #srate at which to interpolate. Default = 4096.
                'use_geometric_units': "yes",   #output quantities in geometric units. Default = "yes"
                'initial_frequency'  : 0.005,    #in Hz if use_geometric_units = 0, else in geometric units
                'interp_uniform_grid': "yes"   #interpolate mode by mode on a uniform grid. Default = "no" (no interpolation)
                #'lal_tetrad_conventions': "no" #for the non-precessing ones in order to get agreement with the modes etc you may have to turn off 
            }
    
teb, hp, hcm, hlm, dyn = EOBRun_module.EOBRunPy(pars)
    
Ah22   = hlm['0'][0] *nu
Phih22 = hlm['0'][1] 
Phih22 = np.unwrap(Phih22)
heb = np.real(Ah22* np.exp(-1j * (Phih22)))
#####

plt.plot(teb,heb,':',label="EOB")
newt = h44.time_ret()/mtot
newh = h44.h.real/mtot
plt.plot(newt-newt[np.argmax(newh)],newh,label='new')
plt.plot(time44-time44[np.argmax(rh44)],rh44,'--',label='old')

plt.legend()
plt.show()
