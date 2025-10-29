from nrtools.initialdata.initialdata import *
from nrtools.evolution.evolution import *
from watpy.wave.wave import mwaves
import os
import numpy as np
import argparse
import matplotlib.pyplot as plt
import glob

## Read args
parser = argparse.ArgumentParser(description='Get Waveform from evolution, e.g. `python get_wvf.py -i id_name -e ev_name -c cluster_name`, either in ARA or DRACO')
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
else:
    print('Error: cluster name unknown')

simpath = os.path.join(basedir,simname)
evo_path = os.path.join(simpath,evo)

## Do evolution(s)
# Initialize ID object
bhns_id = Initial_Data(path=simpath,params=None, id_exe=id_exe)

#ev_folders = [i for i in os.listdir(simpath) if os.path.isdir(os.path.join(simpath,i)) and i.startswith('bam')]
lmax2 = 6
print('\n>> ', evo)
resolution, lmax, flux = ev_folder_parser(evo)
bhns_ev = Evolution(path=os.path.join(simpath,evo), ev_path=ev_path, initial_data=bhns_id, resolution=int(resolution), lmax=int(lmax), lmax2=lmax2,flux=flux) 


ev_output = bhns_ev.ou
id_output = bhns_id.ou
    
mbh, mns, mtot = id_output.get_msun_masses()
_, Momg22 = id_output.get_gw_freqs()

f0 = Momg22 / (2*np.pi) / mtot 

dfiles = [os.path.split(x)[1] for x in glob.glob('{}/{}'.format(ev_output.out_inv_dir,'Rpsi4mode??_r*.l0'))]
wm = mwaves(path = ev_output.out_inv_dir, code = 'bam', filenames = dfiles, 
            mass = mtot, f0 = f0,
            ignore_negative_m=True)

for r in wm.radii:
    for (l,m) in wm.modes:
        hlm = wm.get(l=l, m=m)
        hlm.write_to_txt('h', '.')
#bhns_ev.output_metadata()
