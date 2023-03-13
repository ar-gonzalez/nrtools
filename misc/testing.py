from nrtools.initialdata import *
from nrtools.utils.utils import num_orbits_1pn
import os

# Test with existing ID

outdir = '/beegfs/mo63poy/BHNS_Elliptica/SLy_BH_m2.4_s0--NS_m1.6_s0--d25'

bhns_id = initialdata.Initial_Data(path=outdir,params=None,id_exe=None)

metadata = bhns_id.md

metadic = metadata.id_dic

print(metadic['simname'])

metadata.read_md_sr() # get dic from 2nd highest resolution 

metadic2 = metadata.id_dic2

minf = float(metadic['BH_irreducible_mass_current']) + float(metadic['NS_ADM_mass'])

metadic['BHNS_binding_energy_error'] = str(float(metadic['BHNS_binding_energy'])/minf - float(metadic2['BHNS_binding_energy'])/minf)

metadic['Num_orbits_1PN'] = str(num_orbits_1pn(float(metadic['BH_irreducible_mass_current']),float(metadic['NS_ADM_mass']),float(metadic['BHNS_angular_velocity']),float(metadic['BHNS_ADM_mass'])))

print(metadic)
