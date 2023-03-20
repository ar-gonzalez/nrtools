from nrtools.initialdata.initialdata import *
from nrtools.utils.utils import num_orbits_1pn
import os

# Test with existing ID

outdir = '/home/agonzalez/Documents/PhD/NR/BHNS_evs/evolutions/SLy_BH_m2.8_s0.6--NS_m1.6_s0--d80'
#outdir = '/beegfs/mo63poy/BHNS_Elliptica/SLy_BH_m2.4_s0--NS_m1.6_s0--d25'

bhns_id = Initial_Data(path=outdir,params=None,id_exe='/path/to/elliptica/exe')

metadata = bhns_id.ou

metadic = metadata.id_dic

metadata.read_md_sr() # get dic from 2nd highest resolution 

metadic2 = metadata.id_dic2

minf = float(metadic['BH_irreducible_mass_current']) + float(metadic['NS_ADM_mass'])

metadic['BHNS_binding_energy_error'] = str(float(metadic['BHNS_binding_energy'])/minf - float(metadic2['BHNS_binding_energy'])/minf)

metadic['Num_orbits_1PN'] = str(num_orbits_1pn(float(metadic['BH_irreducible_mass_current']),float(metadic['NS_ADM_mass']),float(metadic['BHNS_angular_velocity']),float(metadic['BHNS_ADM_mass'])))

#bhns_id.check_accuracy()

#bhns_id.check_convergence()

bhns_ev = bhns_id.evolve(ev_path='path/to/bam', resolution=64, lmax=10, lmax2=6, flux='LLF')
bhns_ev.write_bashfile(bashname = 'run_bam.sh', cluster='DRACO')